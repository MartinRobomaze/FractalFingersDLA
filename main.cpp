#include <iostream>
#include <cmath>
#include <random>
#include <fstream>
#include "opencv2/opencv.hpp"
#include "point.h"
#include "convex_hull.h"

#define N_0 (M_PI / (4 * ::sqrt(3)))

using namespace std;

double points_dist(Point p1, Point p2) {
    return sqrt(pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2));
}

int n_points_in_circle(vector<Point> &cluster, Point center, double r) {
    int n_points = 0;

    for (auto p : cluster) {
        if (points_dist(p, center) <= r) {
            n_points++;
        }
    }

    return n_points;
}

pair<double, pair<double, double>> cluster_radius(vector<Point> &cluster) {
    vector<Point> hull;

    convexHull(cluster, hull);

    double maxDistSq = 1;
    pair<double, double> offset_coords;

    for (auto p1 : hull) {
        for (auto p2 : hull) {
            double distSq = pow(p2.x - p1.x, 2) + pow(p2.y - p1.y, 2);

            if (distSq > maxDistSq) {
                offset_coords = pair<double, double>{(p2.x + p1.x) / 2, (p2.y + p1.y) / 2};
                maxDistSq = distSq;
            }
        }
    }

    return pair<double, pair<double, double>>{sqrt(maxDistSq) / 2, offset_coords};
}

int main(int argc, char *argv[]) {
    if (argc != 7) {
        cout << "Usage:" << endl;
        cout << "DLA [R_initial] [num_initial_particles] [num_simulated_particles] [R_diff_circle] [A] [B]" << endl;
        cout << "Number of provided arguments: " << argc << endl;
        return 1;
    }

    double R_initial = stod(argv[1]);
    double num_initial_particles = stod(argv[2]);
    double num_simulated_particles = stod(argv[3]);
    double R_diff_circle = stod(argv[4]);
    double A = stod(argv[5]);
    double B = stod(argv[6]);

    random_device rd;
    mt19937 mt(rd());
    uniform_real_distribution<double> distrib(0.0, 1.0);

    // Initialize the cluster with a single particle at the center
    vector<Point> cluster;
    for (int i = 0; i < num_initial_particles; i++) {
        double angle = 2 * M_PI / num_initial_particles * i;
        double x = R_initial * cos(angle);
        double y = R_initial * sin(angle);

        cluster.push_back(Point{x, y});
    }

    // Set the number of particles and the maximum distance

    double R_diff_circle_sq = pow(R_diff_circle, 2);

    bool wrote_already = false;
    // Generate the particles
    for (int i = 0; i < num_simulated_particles; i++) {
        if (i % 10 == 0 && !wrote_already) {
            cout << i << endl;
            wrote_already = true;
        } else if (i % 10 != 0) {
            wrote_already = false;
        }
        pair<double, pair<double, double>> cluster_R = cluster_radius(cluster);
        // Start the particle at a random distance and angle from the center
        double distance = distrib(mt) * (R_diff_circle) + cluster_R.first;
        double angle = distrib(mt) * 2 * M_PI;
        double x = distance * cos(angle) + cluster_R.second.first;
        double y = distance * sin(angle) + cluster_R.second.second;

        double diff_circle_center_coeff = 1 - R_diff_circle / sqrt(pow(x, 2) + pow(y, 2));
        double x_diff_circle = diff_circle_center_coeff * x;
        double y_diff_circle = diff_circle_center_coeff * y;


        // Random walk the particle until it reaches the cluster
        int edge_reached_count = 0;
        bool alive = true;

        while (alive) {
            double x_prev = x;
            double y_prev = y;
            x += (distrib(mt) - 0.5);
            y += (distrib(mt) - 0.5);

            double sticking_probability =
                    A * ((double)(n_points_in_circle(cluster, Point{x, y}, R_diff_circle)) /
                    (M_PI * pow(R_diff_circle, 2)) - N_0) + B;

            if (sticking_probability < 0.01) {
                sticking_probability = 0.01;
            }

//            cout << n_points_in_circle(cluster, Point{x, y}, R_diff_circle) << " " << sticking_probability << " " << N_0 << endl;

            double diff_circle_dist_sq = pow(x - x_diff_circle, 2) + pow(y - y_diff_circle, 2);
            if (diff_circle_dist_sq > R_diff_circle_sq) {
                edge_reached_count++;
                x = x_prev;
                y = y_prev;
                continue;
            }

            if (edge_reached_count > 2) {
                i--;
                break;
            }

            // Check if the particle is close enough to the cluster.
            for (auto c : cluster) {
                double dx = x - c.x;
                double dy = y - c.y;
                double dist_sq = dx * dx + dy * dy;
                if (dist_sq < 1) {
                    if (distrib(mt) <= sticking_probability) {
                        alive = false;

//                        double dist = sqrt(dist_sq);

                        // Move particle next to other particle.
                        x = c.x + dx;
                        y = c.y + dy;

                        cluster.push_back(Point{x, y});
                    }

                    break;
                }
            }
        }
    }

    pair<double, pair<double, double>> cluster_R = cluster_radius(cluster);

    cout << "Simulation done, saving data" << endl;

    string filename =
            to_string(R_initial) + "_" + to_string(R_diff_circle) + "_" + to_string(A) + "_" +
            to_string(B) + "_" + to_string(num_simulated_particles) + ".csv";

    ofstream file;
    file.open(filename);

    if (!file) {
        cout << "Well, you are fucked..." << endl;
        return -1;
    }

    for (auto c : cluster) {
        file << c.x << "," << c.y << endl;
    }

    file.close();
    return 0;
}
