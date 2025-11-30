Point(1) = {0, 0, 0, 0.001};
Point(2) = {1, 0, 0, 0.001};
Point(3) = {0, 0.025, 0, 0.001};
Line(1)={1,2}; Line(2)={2,3}; Line(3)={3,1};
Line Loop(1) = {1,2,3}; Plane Surface(1) = {1};
