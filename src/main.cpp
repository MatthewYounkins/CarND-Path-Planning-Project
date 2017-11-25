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

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * 0.0174532925199433; } //pi/180       							
double rad2deg(double x) { return x * 57.2957795130823; } //180/pi 

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

double distance(double x1, double y1, double x2, double y2) {
	return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}
int ClosestWaypoint(double x, double y, vector<double> maps_x, vector<double> maps_y) {
	double closestLen = 100000; //large number
	int closestWaypoint = 0;
	for(int i = 0; i < maps_x.size(); i++) {
		double map_x = maps_x[i];
		double map_y = maps_y[i];
		double dist = distance(x,y,map_x,map_y);
		if(dist < closestLen) {
			closestLen = dist;
			closestWaypoint = i;
		}
	}
	return closestWaypoint;
}

int NextWaypoint(double x, double y, double theta, vector<double> maps_x, vector<double> maps_y) {
	int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);
	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];
	double heading = atan2( (map_y-y),(map_x-x) );
	double angle = abs(theta-heading);
	if(angle > 0.785398163) { 
		closestWaypoint++;
	}
	return closestWaypoint;
}

vector<double> getFrenet(double x, double y, double theta,
                         vector<double> maps_x, vector<double> maps_y) {
	int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);
	int prev_wp;
	prev_wp = next_wp-1;
	if(next_wp == 0) 	{
		prev_wp  = maps_x.size()-1;
	}
	double n_x = maps_x[next_wp]-maps_x[prev_wp];
	double n_y = maps_y[next_wp]-maps_y[prev_wp];
	double x_x = x - maps_x[prev_wp];
	double x_y = y - maps_y[prev_wp];
	double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);  	// find the projection of x onto n
	double proj_x = proj_norm*n_x;
	double proj_y = proj_norm*n_y;
	double frenet_d = distance(x_x,x_y,proj_x,proj_y);
	double center_x = 1000-maps_x[prev_wp];
	double center_y = 2000-maps_y[prev_wp];
	double centerToPos = distance(center_x,center_y,x_x,x_y);
	double centerToRef = distance(center_x,center_y,proj_x,proj_y);

	if(centerToPos <= centerToRef) {  
		frenet_d *= -1;
	}

	// calculate s value
	double frenet_s = 0;
	for(int i = 0; i < prev_wp; i++) {
		frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
	}
	frenet_s += distance(0,0,proj_x,proj_y);
	return {frenet_s,frenet_d};
}

// Transforms from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, vector<double> maps_s,
                     vector<double> maps_x, vector<double> maps_y)
{
	int prev_wp = -1;
	while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) )) 	{
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

int main()
{
  uWS::Hub h;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  string map_file_ = "../data/highway_map.csv";		  // Waypoint map to read from
  double max_s = 6945.554;		                      // The max s value before wrapping around the track back to 0

  ifstream in_map_(map_file_.c_str(), ifstream::in);
	std::cout << round(2.3) ;
  string line;
  while (getline(in_map_, line))
	{
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
	double ref_vel = 4;
	time_t startTime = time(0);
	time_t elapsedTime = time(0);
	time_t deltaTime = time(0);
	time_t now = time(0);
	int 	 counter =  0;
	double accel = 0;
	double oldAccel = 0;
	double oldVel = 0;
	double jerk = 0;
	bool   verbose = false;
	
  h.onMessage([&now, &startTime, &elapsedTime, &deltaTime, &accel, &oldAccel,
               &oldVel, &jerk, &counter, 
               &ref_vel, &map_waypoints_x,&map_waypoints_y,
							 &map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy,
							 &lane, &verbose] (uWS::WebSocket<uWS::SERVER> ws, char *data, 														 
							 size_t length, uWS::OpCode opCode) {
	  
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {  // "42" is the answer to everything.  The 4 translates to websocket.  The 2 translates to a websocket event
      auto s = hasData(data);
			if (s != "") {
        auto j = json::parse(s);
        string event = j[0].get<string>();
        if (event == "telemetry") {
          // j[1] is the data JSON object
        	// Our vehicle's Data
        	double carX = j[1]["x"];
         	double carY = j[1]["y"];
         	double carS = j[1]["s"];
         	double carD = j[1]["d"];
         	double carYaw = j[1]["yaw"];
         	double carSpeed = j[1]["speed"];
					double vel = 0.44704*carSpeed; //.44704 = m/s in a MPH
					double desFollowDist = 12;
					double maxAccelRef = 10;	
					double desSpeed;
					double trailModeSpeed;
         	bool 	 trailMode = false;
					bool   tailgating = false;
					bool   lane0IsOK = true;
					bool   lane1IsOK = true;
         	bool   lane2IsOK = true;
					bool   lane3IsOK = false;   // driving on the berm works!
					
					deltaTime = time(0)-now;    //less useful than desired.  This function only gives seconds, we need ms.
					now = time(0);
					elapsedTime = now-startTime;
					counter++;
					accel = 15*(vel-oldVel);  // about 15 iterations per second 
					jerk =  15*(accel-oldAccel);
					oldVel = vel;
					oldAccel = accel;
					
					std::cout << std::fixed;
					std::cout << std::setprecision(2);
			
         	// Previous path data 
         	auto previousPathX = j[1]["previous_path_x"];// Previous paths X values
         	auto previousPathY = j[1]["previous_path_y"];// Previous paths Y values
         	double endPathS = j[1]["end_path_s"];				 // Previous path's end s values 
         	double endPathD = j[1]["end_path_d"];        // Previous path's end d values
         	auto sensor_fusion = j[1]["sensor_fusion"];  // Sensor Data, with all other cars on this side of road.  Head on collisions are avoided by not driving on wrong side of road.
         	int prev_size = previousPathX.size();

					
         	if (prev_size > 0) {
         		carS = endPathS;
         	}

         	for (int i=0; i<sensor_fusion.size(); i++) { // cycle through each car on the road, check each lane to see what's available

						double otherVx = sensor_fusion[i][3];
         		double otherVy = sensor_fusion[i][4];
						double otherS = sensor_fusion[i][5];
						double otherD = sensor_fusion[i][6];
						double otherVel = sqrt(otherVx*otherVx+otherVy*otherVy); 
						double relativeVel = otherVel - vel;				// this would come in handy for some slightly more sophisticated algorithms
						int    otherLane = round(.25*otherD-.5);	
						bool   startThinkingAboutOtherLanes = abs(otherS-carS);
						bool 	 tooClose = abs(otherS-carS) < 30;  	// 22 is tailgating distance, abs(s-carS) is distance between vehicles
						bool   otherCarInFront = otherS-carS > -15;  
						
						if(otherLane == lane && tooClose  ) {
							tailgating = true;
							trailModeSpeed = 2.5*otherVel;  //2.237 MPH in a m/s.  But otherVel / carSpeed seems to be about 2.5? 
							if (verbose) {
								std::cout << "Tailgating:  ";
							}
						}							
						
         		if( otherLane == 0 && tooClose ) { 	// if the car is one lane to my left
							lane0IsOK = false;
         		}
						else if( otherLane == 1 && tooClose) {
							lane1IsOK = false;
						}
         		else if( otherLane == 2 && tooClose) {
         			lane2IsOK = false;
         		}
						
						if(verbose)  { 
							std::cout << i << " " << lane << " " << lane0IsOK << " " << lane1IsOK << " " << lane2IsOK << " "  
						          << carD << " " << carS << " " << ref_vel << " "<< otherS << " " << otherD << " " << otherLane << " " 
											<< relativeVel << " " << otherVel << " " << vel << " " << trailModeSpeed << " " 
											<< carSpeed << " " << tooClose << " " << otherCarInFront << "\n"; 
						}
					}
					if (tailgating) {														//if the previous for loop found a car in the lane, do something
          	if (lane == 0 && lane1IsOK ) {
							lane = 1;
						}
						else if (lane == 1 && lane0IsOK) {
							lane = 0;
						}
						else if (lane == 1 && lane2IsOK) {
							lane = 2;
						}
						else if (lane == 2 && lane1IsOK) {
							lane = 1;
						}
						else {
							trailMode = true;
						}
          }
          

          if (trailMode == false) {
						if(counter < 2) {
								ref_vel = 6;
						}
						else if(counter < 20) {
								ref_vel = 6 + counter*0.15 + counter*counter*.02;
						}
							
						else if(ref_vel < 38) {
							//double coefficient = 2.23694*deltaTime; //this is what it should be... m/s to MPH times deltaTime
							double coefficient = .08;
							ref_vel = coefficient * maxAccelRef + ref_vel; //should be actual velocity but that only gets updated every 15 or so timesteps
						}
						else if (ref_vel < 49.7) {
							ref_vel = ref_vel + 0.1;
						}
						else if (ref_vel > 49.9 ) {
							ref_vel = 49.8;
						}
					}

					if (trailMode == true)  {
						if(verbose) { 
							std::cout << "trailMode: " << trailModeSpeed << " " << carSpeed << " " << ref_vel ; 
							
						}
						if (carSpeed>(trailModeSpeed - .2)) { // slow down to follow the leading car
							//ref_vel = 0.5*(carSpeed + (trailModeSpeed-.2));
							ref_vel = ref_vel - .2;
							std::cout << " Slow to " << ref_vel;
						}/*
						else {  //could mean acceleration in some corner cases
							ref_vel = trailModeSpeed-.1;
							std::cout << "\n Match " << ref_vel << "\n";
						}							*/
         	}
          	
					std::cout << counter << " " << elapsedTime << " " << deltaTime << " " << lane << " " 
										<< carD << " " << carS << " " << trailMode << " " << tailgating << " "
					          << vel << " " << accel << " " << jerk << " " << ref_vel << "\n";

					// create a list of widely spaced (x,y) waypoints, evenly spaced at 30m
					vector<double> ptsx;
					vector<double> ptsy;

					// reference x, y, yaw states
					double ref_x = carX;
					double ref_y = carY;
					double ref_yaw = deg2rad(carYaw);

					if(prev_size < 2)	{
						double prev_carX = carX - cos(carYaw);
						double prev_carY = carY - sin(carYaw);

						ptsx.push_back(prev_carX);
						ptsx.push_back(carX);
						ptsy.push_back(prev_carY);
						ptsy.push_back(carY);
					}
					else {						
						ref_x = previousPathX[prev_size-1];
						ref_y = previousPathY[prev_size-1];	

						double ref_x_prev = previousPathX[prev_size-2];
						double ref_y_prev = previousPathY[prev_size-2];
						ref_yaw = atan2(ref_y-ref_y_prev, ref_x-ref_x_prev);

						ptsx.push_back(ref_x_prev);
						ptsx.push_back(ref_x);

						ptsy.push_back(ref_y_prev);
						ptsy.push_back(ref_y);
					}

					// append points
					vector<double> next_wp0 = getXY(carS+30, (2+4*lane),map_waypoints_s, map_waypoints_x, map_waypoints_y);
					vector<double> next_wp1 = getXY(carS+60, (2+4*lane),map_waypoints_s, map_waypoints_x, map_waypoints_y);
					vector<double> next_wp2 = getXY(carS+90, (2+4*lane),map_waypoints_s, map_waypoints_x, map_waypoints_y);

					ptsx.push_back(next_wp0[0]);
					ptsx.push_back(next_wp1[0]);
					ptsx.push_back(next_wp2[0]);

					ptsy.push_back(next_wp0[1]);
					ptsy.push_back(next_wp1[1]);
					ptsy.push_back(next_wp2[1]);

					for (int i=0; i<ptsx.size(); i++) { // move car reference to zero degrees
						double shift_x = ptsx[i]-ref_x;
						double shift_y = ptsy[i]-ref_y;
						ptsx[i] = (shift_x*cos(0-ref_yaw)-shift_y*sin(0-ref_yaw));
						ptsy[i] = (shift_x*sin(0-ref_yaw)+shift_y*cos(0-ref_yaw));
					}

					tk::spline s;																		// spline curve
					s.set_points(ptsx, ptsy);												// fit spline
					vector<double> next_x_vals;
					vector<double> next_y_vals;	 
					
					for (int i=0; i<previousPathX.size(); i++) {		// add the previuos path for a smooth transition
						next_x_vals.push_back(previousPathX[i]);
						next_y_vals.push_back(previousPathY[i]);
					}

					// keep extending the previous path
					double target_x = 30.0; // m
					double target_y = s(target_x);
					double target_dist = sqrt(target_x*target_x + target_y*target_y);
					double x_add_on = 0;

					// 50 more waypoints
					for (int i=1; i<=50-previousPathX.size();i++)	{
						double N = target_dist/(.02*ref_vel/2.24); // number of intervals
						double x_point = x_add_on+target_x/N;
						double y_point = s(x_point);
						double x_ref = x_point;
						double y_ref = y_point;

						x_add_on = x_point;
						// need to translate x and y back to original coordinates
						x_point = x_ref*cos(ref_yaw) - y_ref*sin(ref_yaw);
						y_point = x_ref*sin(ref_yaw) + y_ref*cos(ref_yaw);
						x_point += ref_x;
						y_point += ref_y;
						next_x_vals.push_back(x_point);
						next_y_vals.push_back(y_point);
					}

				json msgJson;				
				msgJson["next_x"] = next_x_vals;
				msgJson["next_y"] = next_y_vals;
				auto msg = "42[\"control\","+ msgJson.dump()+"]";
				//this_thread::sleep_for(chrono::milliseconds(1000));
				ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
        }
			}
      else {    // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
		}
  });

  // We don't need  this since we're not using HTTP but if it's removed the
  // program  
  // doesn't compile   :-(  
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // could this be done more gracefully?
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
  });

  int port = 4567;  //the best port to listen to. 
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
















































































