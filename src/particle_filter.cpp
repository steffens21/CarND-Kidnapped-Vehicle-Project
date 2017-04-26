/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <random>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  // Set the number of particles. Initialize all particles to first position (based on estimates of
  //   x, y, theta and their uncertainties from GPS) and all weights to 1.
  // Add random Gaussian noise to each particle.
  // NOTE: Consult particle_filter.h for more information about this method (and others in this file).

  cout << " ++ in init" << endl;

  num_particles = 10; // Try more later

  // Random Gaussian noise
  default_random_engine gen;
  normal_distribution<double> N_x(0, std[0]);
  normal_distribution<double> N_y(0, std[1]);
  normal_distribution<double> N_theta(0, std[2]);

  for (int i=0; i<num_particles; i++) {
    Particle p = Particle();

    double n_x = N_x(gen);
    double n_y = N_y(gen);
    double n_theta = N_theta(gen);

    p.id = i;
    p.x = x + n_x;
    p.y = y + n_y;
    p.theta = theta + n_theta;
    p.weight = 1;
    particles.push_back(p);

    cout << "   *** particle " << p.id << " initialized to " << p.x << "\t" << p.y << "\t" << p.theta << endl;
  }

  is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
  // Add measurements to each particle and add random Gaussian noise.
  // NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
  //  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
  //  http://www.cplusplus.com/reference/random/default_random_engine/

  cout << " ++ in prediction" << endl;

  // Random Gaussian noise distributions
  default_random_engine gen;
  normal_distribution<double> N_x(0, std_pos[0]);
  normal_distribution<double> N_y(0, std_pos[1]);
  normal_distribution<double> N_theta(0, std_pos[2]);

  for (int i=0; i < num_particles; i++) {
    Particle p = particles[i];

    // make some noise
    double n_x = N_x(gen);
    double n_y = N_y(gen);
    double n_theta = N_theta(gen);

    if ( yaw_rate > 0.001 ) {
      p.x += velocity / yaw_rate * (sin(p.theta + yaw_rate * delta_t) - sin(p.theta));
      p.y += velocity / yaw_rate * (cos(p.theta) - cos(p.theta + yaw_rate * delta_t));
    } else {
      p.x += cos(p.theta) * velocity * delta_t;
      p.y += sin(p.theta) * velocity * delta_t;
    }
    p.x += n_x;
    p.y += n_y;
    p.theta += yaw_rate * delta_t + n_theta;

    cout << "   *** particle " << p.id << " is now at " << p.x << "\t" << p.y << "\t" << p.theta << endl;
  }
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
  // TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
  //   observed measurement to this particular landmark.
  // NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
  //   implement this method and use it as a helper during the updateWeights phase.

  // predicted: Predicted measurements between one particale and all landmarks within sensor range
  // observations: Actual landmark measruements gathered from lidar

  cout << " ++ in dataAssociation" << endl;

  // Nearest neighbor search
  // This is O( m*n ) where m = size of predicted and n = size of observations.
  for (int j = 0; j < predicted.size(); j++) {
    LandmarkObs pred = predicted[j];
    float min_dist = -1;
    float min_id = -1;
    //cout << "  ** pred lm id " << pred.id << endl;
    for (int k = 0; k < observations.size(); k++) {
      LandmarkObs obs = observations[k];
      float dist_to_lm = dist(obs.x, obs.y, pred.x, pred.y);
      //cout << "obs lm id " << obs.id << endl;
      //cout << "distance " << dist_to_lm << endl;
      if (min_id == -1 or dist_to_lm < min_dist) {
	min_id = obs.id;
	min_dist = dist_to_lm;
      }
    }
    // set the id of predicted landmark to the id of the closet map landmark
    pred.id = min_id;
    //cout << "  ** post pred lm id " << pred.id << endl;
  }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
				   std::vector<LandmarkObs> observations, Map map_landmarks) {
  // TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
  //   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
  // NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
  //   according to the MAP'S coordinate system. You will need to transform between the two systems.
  //   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
  //   The following is a good resource for the theory:
  //   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
  //   and the following is a good resource for the actual equation to implement (look at equation 
  //   3.33. Note that you'll need to switch the minus sign in that equation to a plus to account 
  //   for the fact that the map's y-axis actually points downwards.)
  //   http://planning.cs.uiuc.edu/node99.html

  // - predict measurements to map landmarks within sensor range for each particle
  // - use dataAssociation
  // - update weights
  // - normalize weights?

  // TODO: what to use sensor_range for?
  // TODO: Can I ignore std_landmark since it has been applied to the given observations already?  Maybe use it on the map??

  cout << " ++ in updateWeights" << endl;

  // Remember the sum of all weights for normalization
  float sum_weight = 0.0;
  
  for (int i=0; i < num_particles; i++) {
    Particle p = particles[i];

    // loop over observations
    std::vector<LandmarkObs> observations_trans;
    for (int j = 0; j < observations.size(); ++j) {
      LandmarkObs obs = observations[j];
      LandmarkObs obs_trans;

      obs_trans.id = -1;

      // transform coordinate systems
      obs_trans.x = p.x * cos(p.theta) + p.y * sin(p.theta) + obs.x;
      obs_trans.y = p.x * sin(p.theta) + p.y * cos(p.theta) + obs.y;

      observations_trans.push_back(obs_trans);
    }

    std::vector<LandmarkObs> landmarks_in_range;
    for (int k = 0; k < map_landmarks.landmark_list.size(); ++k) {
      Map::single_landmark_s lm = map_landmarks.landmark_list[k];
      float dist_particle_to_lm = dist(p.x, p.y, lm.x_f, lm.y_f);
      if (dist_particle_to_lm < sensor_range) {
	LandmarkObs in_range_lm;
	in_range_lm.id = lm.id_i;
	in_range_lm.x = lm.x_f;
	in_range_lm.y = lm.y_f;
	landmarks_in_range.push_back(in_range_lm);
      }
    }

    cout << "Nbr. of landmark in range " << landmarks_in_range.size() << endl;
    cout << "Nbr. of observations " << observations_trans.size() << endl;
    //cout << observations_trans << endl;
    //cout << landmarks_in_range << endl;
    dataAssociation(observations_trans, landmarks_in_range);

    float new_weight = 1;
    for (int j = 0; j < observations_trans.size(); ++j) {
      // get coordinates of closets landmark
      LandmarkObs obs = observations_trans[j];
      LandmarkObs closets_lm = landmarks_in_range[obs.id];
      float x_lm = closets_lm.x;
      float y_lm = closets_lm.y;

      // update weights using mult-variate Gaussian distribution
      float normalizer = 1 / (2 * M_PI * std_landmark[0] * std_landmark[1]);
      float first_exp  = pow(obs.x - x_lm, 2)/ pow(std_landmark[0], 2);
      float second_exp = pow(obs.y - y_lm, 2)/ pow(std_landmark[1], 2);
      new_weight *= normalizer * exp(-(first_exp + second_exp));
    }
    p.weight = new_weight;

    sum_weight += new_weight;
  }

  // normalize all weights (maybe not even necessary since discrete_distribution could handle that)
  for (int i=0; i < num_particles; i++) {
    Particle p = particles[i];
    p.weight = p.weight / sum_weight;
  } 
}



void ParticleFilter::resample() {
  // Resample particles with replacement with probability proportional to their weight.
  // NOTE: You may find std::discrete_distribution helpful here.
  //   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

  cout << " ++ in resample" << endl;
  
  random_device rd;
  mt19937 gen(rd());
  discrete_distribution<> d(weights.begin(), weights.end());
  std::vector<Particle> new_particles;
  for(int i=0; i < num_particles; i++) {
    new_particles.push_back( particles[d(gen)] );
  }
  particles = new_particles;
}

void ParticleFilter::write(std::string filename) {
  // You don't need to modify this file.
  std::ofstream dataFile;
  dataFile.open(filename, std::ios::app);
  for (int i = 0; i < num_particles; ++i) {
    dataFile << particles[i].x << " " << particles[i].y << " " << particles[i].theta << "\n";
  }
  dataFile.close();
}
