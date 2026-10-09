// ============================================================================
// gmsh 几何模板 —— 理想圆柱谐振腔 (pillbox)
// ----------------------------------------------------------------------------
// 用法: gmsh -3 cavity.geo -o cavity.msh
// 说明: 半径 R、长度 L 为示意值；Physical 组编号需按 Palace/Elmer 配置调整。
// ============================================================================

R  = 0.30;   // 腔半径 [m]
L  = 0.60;   // 腔长度 [m]
lc = 0.03;   // 网格尺寸 [m]

// 底面圆
Point(1) = {0, 0, 0, lc};
Point(2) = {R, 0, 0, lc};
Point(3) = {0, R, 0, lc};
Point(4) = {-R, 0, 0, lc};
Point(5) = {0, -R, 0, lc};
Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 2};
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// 拉伸成圆柱腔
out[] = Extrude {0, 0, L} { Surface{1}; };
// out[0] = 体, out[1] = 顶面, out[2] = 侧面(柱面)

Physical Volume("vacuum") = {out[0]};
// 注意: 侧面编号取决于 gmsh 版本；请用 gmsh GUI 确认后填入
Physical Surface("wall") = {1, out[1], out[2]};

Mesh.MeshSizeMin = lc;
Mesh.MeshSizeMax = lc*2;
Mesh 3;
