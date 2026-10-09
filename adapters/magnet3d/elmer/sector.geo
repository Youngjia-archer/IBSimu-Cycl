// ============================================================================
// gmsh 几何模板 —— 回旋加速器扇形磁铁（示意）
// ----------------------------------------------------------------------------
// 用法: gmsh -3 sector.geo -o sector.msh
// 说明: 这是**占位/示意**几何。实际使用时请按你的磁铁参数（扇形角、磁极半径、
//       气隙高度、回旋角 spiral angle、垫补 shim 等）修改。
//       也支持从 CAD 导入: Merge "yoke.step";
// ============================================================================

// --- 参数 -------------------------------------------------------------------
R_in   = 1.80;   // 内半径 [m]
R_out  = 4.80;   // 外半径 [m]
Sector = 45.0;   // 扇形张角 [deg]  (8 扇形 -> 45)
Gap    = 0.06;   // 上下磁极气隙半高 [m]
PoleH  = 0.40;   // 磁极高度 [m]
YokeH  = 0.40;   // 铁轭厚度 [m]
h      = 0.15;   // 网格尺寸 [m]

// --- 铁矿物体（简化：径向扇块） ---------------------------------------------
// 这里仅示意一个扇块；真实几何需包含磁极、磁轭回路、垫补等
r1 = R_in;  r2 = R_out;
Point(1) = { r1, 0, -Gap, h };
Point(2) = { r2, 0, -Gap, h };
Point(3) = { r1, 0, -Gap-PoleH, h };
Point(4) = { r2, 0, -Gap-PoleH, h };
Point(5) = { r1, 0,  Gap, h };
Point(6) = { r2, 0,  Gap, h };
Point(7) = { r1, 0,  Gap+PoleH, h };
Point(8) = { r2, 0,  Gap+PoleH, h };

Line(1) = {1,2}; Line(2) = {2,4}; Line(3) = {4,3}; Line(4) = {3,1};
Line(5) = {5,6}; Line(6) = {6,8}; Line(7) = {8,7}; Line(8) = {7,5};

Curve Loop(1) = {1,2,3,4};  Plane Surface(1) = {1};   // 下磁极截面
Curve Loop(2) = {5,6,7,8};  Plane Surface(2) = {2};   // 上磁极截面
Curve Loop(3) = {1,-5,8,-4,3,-7,6,-2}; Plane Surface(3) = {3};  // 侧壁(示意)

Surface Loop(1) = {1,2,3};
Volume(1) = {1};                 // 铁矿(示意)

// --- 空气域（包围） ---------------------------------------------------------
Box_air = 6.0;
Point(101) = { -Box_air, -Box_air, -Box_air, h*2 };
Point(102) = {  Box_air, -Box_air, -Box_air, h*2 };
Point(103) = {  Box_air,  Box_air, -Box_air, h*2 };
Point(104) = { -Box_air,  Box_air, -Box_air, h*2 };
Point(105) = { -Box_air, -Box_air,  Box_air, h*2 };
Point(106) = {  Box_air, -Box_air,  Box_air, h*2 };
Point(107) = {  Box_air,  Box_air,  Box_air, h*2 };
Point(108) = { -Box_air,  Box_air,  Box_air, h*2 };

Line(101) = {101,102}; Line(102) = {102,103}; Line(103) = {103,104}; Line(104) = {104,101};
Line(105) = {105,106}; Line(106) = {106,107}; Line(107) = {107,108}; Line(108) = {108,105};
Line(109) = {101,105}; Line(110) = {102,106}; Line(111) = {103,107}; Line(112) = {104,108};

Curve Loop(101) = {101,102,103,104}; Plane Surface(101) = {101};
Curve Loop(102) = {105,106,107,108}; Plane Surface(102) = {102};
Curve Loop(103) = {101,110,-105,-109}; Plane Surface(103) = {103};
Curve Loop(104) = {102,111,-106,-110}; Plane Surface(104) = {104};
Curve Loop(105) = {103,112,-107,-111}; Plane Surface(105) = {105};
Curve Loop(106) = {104,109,-108,-112}; Plane Surface(106) = {106};

Surface Loop(101) = {101,102,103,104,105,106};
Volume(2) = {101};   // 空气

// --- 物理组（对应 .sif 中的 Target Bodies / Target Boundaries） --------------
Physical Volume("air")  = {2};
Physical Volume("iron") = {1};
Physical Surface("farfield") = {101,102,103,104,105,106};

Mesh.MeshSizeMin = h;
Mesh.MeshSizeMax = h*2;
Mesh 3;
