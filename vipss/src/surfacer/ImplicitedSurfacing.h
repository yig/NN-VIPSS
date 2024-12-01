#ifndef IMPLICITEDSURFACING_H
#define IMPLICITEDSURFACING_H

#include "Polygonizer.h"
#include "../readers.h"

typedef unsigned int uint;

class Surfacer{

public:

   R3Pt							s_ptMin, s_ptMax;
   Array<R3Pt>					s_aptSurface;
   Array<R3Vec>					s_avecSurface;
   Array<R3Pt_i>				s_afaceSurface;


   std::vector<double>all_v;
   std::vector<uint>all_fv;

   R3Pt st;
   R3Pt center;
   double dSize;
   int iBound;


   Surfacer(){}

   void CalSurfacingPara(std::vector<double>&Vs, int nvoxels);

   double Surfacing_Implicit(std::vector<double>&Vs, int n_voxels, bool ischeckall,
                  double (*function)(const R3Pt &in_pt));



   void WriteSurface(std::string fname);
   void WriteSurface(std::vector<double> &v, std::vector<uint>&fv);
   void WriteSurface(std::vector<double> **v, std::vector<uint>**fv);

   void ClearBuffer();

   void ClearSingleComponentBuffer();

private:
   void GetCurSurface(std::vector<double> &v, std::vector<uint>&fv);
   void InsertToCurSurface(std::vector<double>&v, std::vector<uint>&fv);





};


















#endif // IMPLICITEDSURFACING_H
