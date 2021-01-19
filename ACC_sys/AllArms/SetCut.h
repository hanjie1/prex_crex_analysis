#include "CollimatorL.C"
#include "CollimatorR.C"
#include "UpPlane.C"

TString xcut = "x_fp_tr!=-333.";
TString colcutL = "CollimatorL(x_col_tr,y_col_tr)";
TString colcutR = "CollimatorR(x_col_tr,y_col_tr)";

TCut XCUT = Form("%s",xcut.Data());//&&x_vdc_tr>0.0";
TCut colCutL = Form("%s",colcutL.Data());
TCut colCutR = Form("%s",colcutR.Data());
TCut UPCutL = "UpPlane(x_zup1,y_zup1,x_zup2,y_zup2,1)";
TCut UPCutR = "UpPlane(x_zup1,y_zup1,x_zup2,y_zup2,0)";

double xmin[9] = { 0.082, 0.0862, 0.0989, 0.1084, 0.1195, 0.1236, 0.1326, 0.1958, 0.2149 };
double xmax[9] = { 0.1256, 0.1318, 0.1505, 0.1643, 0.1798, 0.1853, 0.1975, 0.2777, 0.301 };
double ymin[9] = {-0.0447, -0.0469, -0.053, -0.0571, -0.0612, -0.0624, -0.065, -0.0786, -0.0818 };
double ymax[9] = { 0.0445, 0.0466, 0.0529, 0.057, 0.0611, 0.0624, 0.065, 0.0785, 0.0818 };

double pinch = 0.0023;

TCut downcutL1 = Form("x_zdown1>%f && x_zdown1<%f && y_zdown1>%f && y_zdown1<%f",xmin[0],xmax[0],ymin[0],ymax[0]);
TCut downcutL2 = Form("x_zdown2>%f && x_zdown2<%f && y_zdown2>%f && y_zdown2<%f",xmin[1],xmax[1],ymin[1],ymax[1]);
TCut downcutL3 = Form("x_zdown3-%f>%f && x_zdown3-%f<%f && y_zdown3>%f && y_zdown3<%f",pinch,xmin[2],pinch,xmax[2],ymin[2],ymax[2]);
TCut downcutL4 = Form("x_zdown4-%f>%f && x_zdown4-%f<%f && y_zdown4>%f && y_zdown4<%f",pinch,xmin[3],pinch,xmax[3],ymin[3],ymax[3]);
TCut downcutL5 = Form("x_zdown5>%f && x_zdown5<%f && y_zdown5>%f && y_zdown5<%f",xmin[4],xmax[4],ymin[4],ymax[4]);
TCut downcutL6 = Form("x_zdown6>%f && x_zdown6<%f && y_zdown6>%f && y_zdown6<%f",xmin[5],xmax[5],ymin[5],ymax[5]);
TCut downcutL7 = Form("x_zdown7>%f && x_zdown7<%f && y_zdown7>%f && y_zdown7<%f",xmin[6],xmax[6],ymin[6],ymax[6]);
TCut downcutL8 = Form("x_zdown8>%f && x_zdown8<%f && y_zdown8>%f && y_zdown8<%f",xmin[7],xmax[7],ymin[7],ymax[7]);
TCut downcutL9 = Form("x_zdown9>%f && x_zdown9<%f && y_zdown9>%f && y_zdown9<%f",xmin[8],xmax[8],ymin[8],ymax[8]);

TCut DownCutL =downcutL1+downcutL2+downcutL3+downcutL4+downcutL5+downcutL6+downcutL7+downcutL8+downcutL9;


TCut downcutR1 = Form("-x_zdown1>%f && -x_zdown1<%f && -y_zdown1>%f && -y_zdown1<%f",xmin[0],xmax[0],ymin[0],ymax[0]);
TCut downcutR2 = Form("-x_zdown2>%f && -x_zdown2<%f && -y_zdown2>%f && -y_zdown2<%f",xmin[1],xmax[1],ymin[1],ymax[1]);
TCut downcutR3 = Form("-x_zdown3>%f && -x_zdown3<%f && -y_zdown3>%f && -y_zdown3<%f",xmin[2],xmax[2],ymin[2],ymax[2]);
TCut downcutR4 = Form("-x_zdown4>%f && -x_zdown4<%f && -y_zdown4>%f && -y_zdown4<%f",xmin[3],xmax[3],ymin[3],ymax[3]);
TCut downcutR5 = Form("-x_zdown5>%f && -x_zdown5<%f && -y_zdown5>%f && -y_zdown5<%f",xmin[4],xmax[4],ymin[4],ymax[4]);
TCut downcutR6 = Form("-x_zdown6>%f && -x_zdown6<%f && -y_zdown6>%f && -y_zdown6<%f",xmin[5],xmax[5],ymin[5],ymax[5]);
TCut downcutR7 = Form("-x_zdown7>%f && -x_zdown7<%f && -y_zdown7>%f && -y_zdown7<%f",xmin[6],xmax[6],ymin[6],ymax[6]);
TCut downcutR8 = Form("-x_zdown8>%f && -x_zdown8<%f && -y_zdown8>%f && -y_zdown8<%f",xmin[7],xmax[7],ymin[7],ymax[7]);
TCut downcutR9 = Form("-x_zdown9>%f && -x_zdown9<%f && -y_zdown9>%f && -y_zdown9<%f",xmin[8],xmax[8],ymin[8],ymax[8]);

TCut DownCutR =downcutR1+downcutR2+downcutR3+downcutR4+downcutR5+downcutR6+downcutR7+downcutR8+downcutR9;


// lead cut
TString ispb = "ev.nuclA==208";
TCut isPb = Form("%s",ispb.Data());
