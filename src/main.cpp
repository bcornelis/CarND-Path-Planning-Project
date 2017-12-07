#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include "spline.h"

using namespace std;

double car_acceptable_boost=.224;
double const legal_speed_limit = 49.5;
double speed_lookahead = 30;
double behaviour_cost_lookahead_distance=speed_lookahead+5;
int prediction_steps = 30;

double approx_speed = car_acceptable_boost;

// for convenience
using json = nlohmann::json;

class Car {

public:
	double x,y,yaw,prev_speed;
	double s,d;
	int _lane;
	double _next_speed;

	Car(const json& input, const int lane) {
		x = input["x"];
		y = input["y"];
		s = input["s"];
		d = input["d"];
		yaw = input["yaw"];
		prev_speed = input["speed"];
		_lane = lane;
	}

	void setNextSpeed(const double nspeed) {
		_next_speed = nspeed;
	}
};

struct Road {
	json previous_path_x, previous_path_y;
	double end_path_s, end_path_d;

	Road(const json& input) {
		 previous_path_x = input["previous_path_x"];
		 previous_path_y = input["previous_path_y"];
		 // Previous path's end s and d values
		 end_path_s = input["end_path_s"];
		 end_path_d = input["end_path_d"];
	}
};

struct Prediction {
	int s;
	double d;
	double speed;
};



// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

void apply_trajectory_KL(map<int, vector<Prediction>>& predictions, Car& car, const Road& road) {

	// set to maximum value. All functions are allowed to lower, but not to increase!
	car.setNextSpeed(10000);

	// iterate over all cars
	for(int carIdx=0; carIdx<predictions.size(); carIdx++) {
		float other_d = predictions[carIdx][0].d;
		// does the car collides with our car?
		if(other_d < (2+4*car._lane+2) && other_d > (2+4*car._lane-2)) {
			//double vx = sensor_fusion[carIdx][3];
			//double vy = sensor_fusion[carIdx][4];
			//double check_speed = sqrt(vx*vx+vy*vy);
			double check_speed = predictions[carIdx][0].speed;
			double check_car_s = predictions[carIdx][0].s;
			double distance_to_car = check_car_s - car.s;

			// where will the car be in the future?
			if( (check_car_s > car.s) && distance_to_car < speed_lookahead) {
				// are we driving approximately the same speed as the other car?
				car.setNextSpeed(check_speed/0.44704);
			}
		}
	}
}

bool is_car_in_the_way_for_lane_change(map<int, vector<Prediction>>& predictions, Car& car, const int& lane) {
	double my_car_behind_s = car.s - 2.;
	double my_car_front_s = car.s + 20.;

	// is there a car int he lane?
	for(int carIdx=0; carIdx<predictions.size(); carIdx++) {
		float other_d = predictions[carIdx][0].d;
		// does the car collides with our car?
		if(other_d < (2+4*lane+2) && other_d > (2+4*lane-2)) {

			for( int posCtr=0; posCtr<10; posCtr++) {
				double check_car_s = predictions[carIdx][posCtr].s;

				if( check_car_s > my_car_behind_s && check_car_s < my_car_front_s) {
					return true;
				}
			}
		}
	}
	return false;
}

void apply_trajectory_LCL(map<int, vector<Prediction>>& predictions, Car& car, const Road& road, int& lane) {

	if( is_car_in_the_way_for_lane_change(predictions, car, lane-1) ) {
		apply_trajectory_KL(predictions, car, road);
	} else {
		lane -=1;
	}
}

void apply_trajectory_LCR(map<int, vector<Prediction>>& predictions, Car& car, const Road& road, int& lane) {

	if( is_car_in_the_way_for_lane_change(predictions, car, lane+1) ) {
		apply_trajectory_KL(predictions, car, road);
	} else {
		lane +=1;
	}
}

double distance(double x1, double y1, double x2, double y2)
{
	return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}
int ClosestWaypoint(double x, double y, const vector<double> &maps_x, const vector<double> &maps_y)
{

	double closestLen = 100000; //large number
	int closestWaypoint = 0;

	for(int i = 0; i < maps_x.size(); i++)
	{
		double map_x = maps_x[i];
		double map_y = maps_y[i];
		double dist = distance(x,y,map_x,map_y);
		if(dist < closestLen)
		{
			closestLen = dist;
			closestWaypoint = i;
		}

	}

	return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{

	int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];

	double heading = atan2( (map_y-y),(map_x-x) );

	double angle = abs(theta-heading);

	if(angle > pi()/4)
	{
		closestWaypoint++;
	}

	return closestWaypoint;

}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

	int prev_wp;
	prev_wp = next_wp-1;
	if(next_wp == 0)
	{
		prev_wp  = maps_x.size()-1;
	}

	double n_x = maps_x[next_wp]-maps_x[prev_wp];
	double n_y = maps_y[next_wp]-maps_y[prev_wp];
	double x_x = x - maps_x[prev_wp];
	double x_y = y - maps_y[prev_wp];

	// find the projection of x onto n
	double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
	double proj_x = proj_norm*n_x;
	double proj_y = proj_norm*n_y;

	double frenet_d = distance(x_x,x_y,proj_x,proj_y);

	//see if d value is positive or negative by comparing it to a center point

	double center_x = 1000-maps_x[prev_wp];
	double center_y = 2000-maps_y[prev_wp];
	double centerToPos = distance(center_x,center_y,x_x,x_y);
	double centerToRef = distance(center_x,center_y,proj_x,proj_y);

	if(centerToPos <= centerToRef)
	{
		frenet_d *= -1;
	}

	// calculate s value
	double frenet_s = 0;
	for(int i = 0; i < prev_wp; i++)
	{
		frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
	}

	frenet_s += distance(0,0,proj_x,proj_y);

	return {frenet_s,frenet_d};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, const vector<double> &maps_s, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int prev_wp = -1;

	while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
	{
		prev_wp++;
	}

	int wp2 = (prev_wp+1)%maps_x.size();

	double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
	// the x,y,s along the segment
	double seg_s = (s-maps_s[prev_wp]);

	double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
	double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

	double perp_heading = heading-pi()/2;

	double x = seg_x + d*cos(perp_heading);
	double y = seg_y + d*sin(perp_heading);

	return {x,y};

}

map<int, vector<Prediction>> predict(const json& sensor_fusion, int prev_size) {

	std::map<int, vector<Prediction>> retval;

	for(int carIdx=0; carIdx<sensor_fusion.size(); carIdx++) {

		// get the cars properties
		double vx = sensor_fusion[carIdx][3];
		double vy = sensor_fusion[carIdx][4];
		double check_speed = sqrt(vx*vx+vy*vy);
		double check_car_s = sensor_fusion[carIdx][5];
		float other_d = sensor_fusion[carIdx][6];

		vector<Prediction> positions;

		for(int prd_step=prev_size; prd_step<50; prd_step++) {
			double s = check_car_s + (prd_step * 0.02 * check_speed);

			Prediction newPrediction;
			newPrediction.s = s;
			newPrediction.d = other_d;
			newPrediction.speed = check_speed;
			positions.push_back(newPrediction);

		}

		// and add to the map
		retval[carIdx] = positions;
	}

	return retval;
}

/**
 * The idea is: the shorter the smallest distance between the ego car and any other car in the lane is
 * the cost. The shorter the distance, the higher the cost
 */
double cost_for_lane(const Car& myself, const int lane, map<int, vector<Prediction>>& predictions) {

	double shortest_car_distance_front = numeric_limits<double>::infinity();
	double shortest_car_speed_front = numeric_limits<double>::infinity();

	// based on speed
	for(int carIdx=0; carIdx<predictions.size(); carIdx++) {
		// we assume the lane is constant in the predictions
		float other_d = predictions[carIdx][0].d;

		// is the car in the required lane?
		if(other_d < (2+4*lane+2) && other_d > (2+4*lane-2)) {
			// yes, find the first car in front, and check distance and speed
			double check_speed = predictions[carIdx][0].speed;
			double check_car_s = predictions[carIdx][0].s;
			double distance_to_car = check_car_s - myself.s;

			// where will the car be in the future?
			if( (check_car_s > myself.s) && distance_to_car < behaviour_cost_lookahead_distance && distance_to_car < shortest_car_distance_front) {
				// newest shortest in front
				shortest_car_distance_front = distance_to_car;
				shortest_car_speed_front = check_speed;
			}
		}
	}

	// if we didn't find a result
	if( shortest_car_distance_front > 10000000.) {
		// cost is optimal
		return 0;
	} else {
		double cost_for_distance = 1. - (shortest_car_distance_front/behaviour_cost_lookahead_distance);
		double cost_for_speed = 1. - shortest_car_speed_front / 60.;
		//return  cost_for_distance * cost_for_speed;
		return cost_for_speed;
	}
}

// TODO: car must be const!
void moveCar(Car& car, const Road& road, const vector<double>& map_waypoints_x, const vector<double>& map_waypoints_y,
		const vector<double>& map_waypoints_s, const vector<double>& map_waypoints_dx, const vector<double>& map_waypoints_dy,
		vector<double>& next_x_vals, vector<double>& next_y_vals) {

	// some easy references
	int prev_size = road.previous_path_x.size();

	// create a list of from which to generate a spline
	vector<double> ptsx;
	vector<double> ptsy;

	//the current cars coordinates
	double ref_x = car.x;
	double ref_y = car.y;
	double ref_yaw = deg2rad(car.yaw);

	// if the list of previous points is (almost) empty, use the car as the starting reference
	// we store those 2 points to have the proper starting point for our spline, and the proper initial angle
	if( prev_size < 2 ) {
		// use the current car's position as the reference
		double prev_car_x = car.x - cos(car.yaw);
		double prev_car_y = car.y - sin(car.yaw);

		ptsx.push_back(prev_car_x);
		ptsx.push_back(car.x);
		ptsy.push_back(prev_car_y);
		ptsy.push_back(car.y);
	} else {
		// there's enough points left in the array of previously generated points
		// select the last one as the current position
		ref_x = road.previous_path_x[prev_size-1];
		ref_y = road.previous_path_y[prev_size-1];

		// and the one before latest as the position one t before
		double ref_x_prev = road.previous_path_x[prev_size-2];
		double ref_y_prev = road.previous_path_y[prev_size-2];
		// and calculate the 'current' yaw based on the two last points
		ref_yaw = atan2(ref_y-ref_y_prev, ref_x-ref_x_prev);

		// store the two latest points in our array
		ptsx.push_back(ref_x_prev);
		ptsx.push_back(ref_x);
		ptsy.push_back(ref_y_prev);
		ptsy.push_back(ref_y);
	}

	// we now have the two latest points in our array. To generate next points the steps are:
	// 1. create 3 points 'far' away (they keep lane changes,speed etc into account. We have one point 30meters ahead, one 60 meters ahead and one 90 meters ahead
	vector<double> next_mp0 = getXY(car.s+30, (2+4*car._lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);
	vector<double> next_mp1 = getXY(car.s+50, (2+4*car._lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);
	vector<double> next_mp2 = getXY(car.s+70, (2+4*car._lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);
	vector<double> next_mp3 = getXY(car.s+90, (2+4*car._lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);

	// store them into our array of path generation points
	ptsx.push_back(next_mp0[0]);
	ptsx.push_back(next_mp1[0]);
	ptsx.push_back(next_mp2[0]);
	ptsx.push_back(next_mp3[0]);

	ptsy.push_back(next_mp0[1]);
	ptsy.push_back(next_mp1[1]);
	ptsy.push_back(next_mp2[1]);
	ptsy.push_back(next_mp3[1]);

	// 2. generate a spline through those points. As the first two points contain the proper angle, it will properly
	// 'continue' on the previously generated splines from which the points are already on the previous_path
	// first, to make the math easier afterwards, a transformation is performed to shift the points into the cars coordinate system
	for(int i=0; i<ptsx.size(); i++) {
		// precalculate how much we need to shift (=difference between the current point and the cars location)
		// This also implies we'll have to perform the inverse transformation later on to have the points in the
		// original x/y coordinate system
		double shift_x = ptsx[i]-ref_x;
		double shift_y = ptsy[i]-ref_y;

		// rotation + translation
		ptsx[i] = (shift_x * cos(0-ref_yaw)-shift_y*sin(0-ref_yaw));
		ptsy[i] = (shift_x * sin(0-ref_yaw)+shift_y*cos(0-ref_yaw));
	}

	// generate the spline through our points
	tk::spline s;
	s.set_points(ptsx, ptsy);
	// precalculate some of the parameters to use later on for the extra point generation
	double target_x = 30; // in 30meters...
	double target_y = s(target_x); // what's the spline's y-value for the specified x-value
	double target_dist = sqrt((target_x*target_x)+(target_y*target_y));

	// 3. fill the array with points
	// fill the new list of points with all left-over points of the previously generated points
	for(int i=0; i<road.previous_path_x.size(); i++) {
		next_x_vals.push_back(road.previous_path_x[i]);
		next_y_vals.push_back(road.previous_path_y[i]);
	}

	// and now, generate more points to fill the array
	double x_add_on = 0;
	for(int i=1; i<=50-prev_size; i++) {
		// should we go faster?
		if(car._next_speed > approx_speed) {
			approx_speed += car_acceptable_boost;
		} else {
			approx_speed -= car_acceptable_boost;
		}

		// some final adjustements
		if( approx_speed > legal_speed_limit) {
			approx_speed = legal_speed_limit;
		}

		//std::cout << "OrigSpeed: " << car.prev_speed << " ;speed_loop: " << approx_speed << std::endl;

		double N = target_dist/(.02 * approx_speed / 2.24);
		double x_point = x_add_on + (target_x) / N; // the next x-point
		double y_point = s(x_point); // the next y-point (provided by the spline)
		x_add_on = x_point;

		// remember to transfer the points back from car coordinates into worls coordinates
		double x_ref = x_point;
		double y_ref = y_point;

		x_point = (x_ref * cos(ref_yaw) - y_ref * sin(ref_yaw));
		y_point = (x_ref * sin(ref_yaw) + y_ref * cos(ref_yaw));

		x_point += ref_x;
		y_point += ref_y;

		next_x_vals.push_back(x_point);
		next_y_vals.push_back(y_point);
	}
}

int main() {
  uWS::Hub h;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // Waypoint map to read from
  string map_file_ = "./data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;

  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;
  while (getline(in_map_, line)) {
  	istringstream iss(line);
  	double x;
  	double y;
  	float s;
  	float d_x;
  	float d_y;
  	iss >> x;
  	iss >> y;
  	iss >> s;
  	iss >> d_x;
  	iss >> d_y;
  	map_waypoints_x.push_back(x);
  	map_waypoints_y.push_back(y);
  	map_waypoints_s.push_back(s);
  	map_waypoints_dx.push_back(d_x);
  	map_waypoints_dy.push_back(d_y);
  }

  int lane = 1;

  h.onMessage([&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy,&lane](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);
        
        string event = j[0].get<string>();
        
        if (event == "telemetry") {
          // j[1] is the data JSON object
          
        	// Main car's localization Data
          	/*double car_x = j[1]["x"];
          	double car_y = j[1]["y"];
          	double car_s = j[1]["s"];
          	double car_d = j[1]["d"];
          	double car_yaw = j[1]["yaw"];
          	double car_speed = j[1]["speed"];*/


          	// Previous path data given to the Planner
          	auto previous_path_x = j[1]["previous_path_x"];
          	auto previous_path_y = j[1]["previous_path_y"];
          	// Previous path's end s and d values 
          	double end_path_s = j[1]["end_path_s"];
          	double end_path_d = j[1]["end_path_d"];

          	// Sensor Fusion Data, a list of all other cars on the same side of the road.
          	auto sensor_fusion = j[1]["sensor_fusion"];

          	Car myself = Car(j[1], lane);
          	Road road = Road(j[1]);

          	if( previous_path_x.size() > 0) {
          		myself.s = end_path_s;
          	}

          	json msgJson;

          	vector<double> next_x_vals;
          	vector<double> next_y_vals;


          	// definitions

          	// 1. prediction
          	map<int, vector<Prediction>> predictions = predict(sensor_fusion, previous_path_x.size());

          	// 2. behaviour planning
          	string next_state = "KL";
          	double cost_for_KL = cost_for_lane(myself, lane, predictions);
          	double cost_for_LCL = (lane == 0) ? 5. : cost_for_lane(myself, lane -1, predictions);
          	double cost_for_LCR = (lane == 2) ? 5. : cost_for_lane(myself, lane + 1, predictions);

            if(lane != 2 && cost_for_LCR < cost_for_KL && cost_for_LCR < cost_for_LCL ) {
          		next_state = "LCR";
          	} else if( lane != 0 && cost_for_LCL < cost_for_KL && cost_for_LCL < cost_for_LCR ) {
          		next_state = "LCL";
          	}

          	// 3. trajectory generation
          	if( next_state == "KL" ) apply_trajectory_KL(predictions, myself, road);
          	if( next_state == "LCL") apply_trajectory_LCL(predictions, myself, road, lane);
          	if( next_state == "LCR") apply_trajectory_LCR(predictions, myself, road, lane);
          	moveCar(myself, road, map_waypoints_x, map_waypoints_y, map_waypoints_s, map_waypoints_dx, map_waypoints_dy, next_x_vals, next_y_vals);

          	// TODO: define a path made up of (x,y) points that the car will visit sequentially every .02 seconds
          	msgJson["next_x"] = next_x_vals;
          	msgJson["next_y"] = next_y_vals;

          	auto msg = "42[\"control\","+ msgJson.dump()+"]";

          	//this_thread::sleep_for(chrono::milliseconds(1000));
          	ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);

        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
 /* h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });*/

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
