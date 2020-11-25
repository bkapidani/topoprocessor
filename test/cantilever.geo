// Gmsh project created on Mon Nov 23 23:37:26 2020
Merge "test.brep";
// For the progressive distribution of element size on some edges. HoWill this is for you!
Bias=1.2;
Transfinite Line{1,2,7,9,15,20,21,24} = 3;
Transfinite Line{3,4,11,12,16,18,25,26} = 10;
// Progressive distribution. Some lines are not oriented in the same way. So the invers of Bias is used.
Transfinite Line{22,23} = 10 Using Progression 1/Bias;
Transfinite Line{27,28} = 10 Using Progression Bias;
Transfinite Line{13,14,17,19} = 10;
Transfinite Line{5,6,8,10} = 4;
Transfinite Surface "*";
Recombine Surface "*";
Transfinite Volume "*";
