//
// Created by robomaze on 27.1.2023.
//
#include <vector>
#include "point.h"

using namespace std;

#ifndef DLA_CPP_CONVEX_HULL_H
#define DLA_CPP_CONVEX_HULL_H

Point nextToTop(vector<Point> &S);
void swap(Point &p1, Point &p2);
double distSq(Point p1, Point p2);
int orientation(Point p, Point q, Point r);
int compare(const void *vp1, const void *vp2);
void convexHull(vector<Point> points, vector<Point> &hull);

#endif //DLA_CPP_CONVEX_HULL_H
