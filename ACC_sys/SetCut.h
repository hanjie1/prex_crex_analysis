#include "CollimatorL.C"
#include "UpPlane.C"

TString xcut = "x_fp_tr!=-333.";
TString colcut = "CollimatorL(x_col_tr,y_col_tr)";

TCut XCUT = Form("%s",xcut.Data());//&&x_vdc_tr>0.0";
TCut colCut = Form("%s",colcut.Data());
TCut UPCut = "UpPlane(x_zup1,y_zup1,x_zup2,y_zup2,1)";

double xmin[9] = { 0.082, 0.0862, 0.0989, 0.1084, 0.1195, 0.1236, 0.1326, 0.1958, 0.2149 };
double xmax[9] = { 0.1256, 0.1318, 0.1505, 0.1643, 0.1798, 0.1853, 0.1975, 0.2777, 0.301 };
double ymin[9] = {-0.0447, -0.0469, -0.053, -0.0571, -0.0612, -0.0624, -0.065, -0.0786, -0.0818 };
double ymax[9] = { 0.0445, 0.0466, 0.0529, 0.057, 0.0611, 0.0624, 0.065, 0.0785, 0.0818 };

double pinch = 0.0023;

TCut downcut1 = Form("x_zdown1>%f && x_zdown1<%f && y_zdown1>%f && y_zdown1<%f",xmin[0],xmax[0],ymin[0],ymax[0]);
TCut downcut2 = Form("x_zdown2>%f && x_zdown2<%f && y_zdown2>%f && y_zdown2<%f",xmin[1],xmax[1],ymin[1],ymax[1]);
TCut downcut3 = Form("x_zdown3-%f>%f && x_zdown3-%f<%f && y_zdown3>%f && y_zdown3<%f",pinch,xmin[2],pinch,xmax[2],ymin[2],ymax[2]);
TCut downcut4 = Form("x_zdown4-%f>%f && x_zdown4-%f<%f && y_zdown4>%f && y_zdown4<%f",pinch,xmin[3],pinch,xmax[3],ymin[3],ymax[3]);
TCut downcut5 = Form("x_zdown5>%f && x_zdown5<%f && y_zdown5>%f && y_zdown5<%f",xmin[4],xmax[4],ymin[4],ymax[4]);
TCut downcut6 = Form("x_zdown6>%f && x_zdown6<%f && y_zdown6>%f && y_zdown6<%f",xmin[5],xmax[5],ymin[5],ymax[5]);
TCut downcut7 = Form("x_zdown7>%f && x_zdown7<%f && y_zdown7>%f && y_zdown7<%f",xmin[6],xmax[6],ymin[6],ymax[6]);
TCut downcut8 = Form("x_zdown8>%f && x_zdown8<%f && y_zdown8>%f && y_zdown8<%f",xmin[7],xmax[7],ymin[7],ymax[7]);
TCut downcut9 = Form("x_zdown9>%f && x_zdown9<%f && y_zdown9>%f && y_zdown9<%f",xmin[8],xmax[8],ymin[8],ymax[8]);

TCut DownCut =downcut1+downcut2+downcut3+downcut4+downcut5+downcut6+downcut7+downcut8+downcut9;

// lead cut
TString ispb = "ev.nuclA==208";
TCut isPb = Form("%s",ispb.Data());
