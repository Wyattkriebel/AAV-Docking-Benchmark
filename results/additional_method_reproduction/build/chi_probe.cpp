#include <iostream>
#include <iomanip>
#include <stdexcept>
#include "parse_pdbqt.h"
#include "model.h"
int main(int argc,char** argv) {
 if(argc!=2) return 2;
 model m=parse_ligand_pdbqt(path(argv[1]));
 vecv coordinates=m.get_flexible_coords();
 std::vector<parsed_atom> input=liginfo_return();
 std::vector<std::vector<size_t*> > glyco=glycan_info_func();
 size_t supported=0;
 for(size_t j=0;j<glyco.size();++j) {
  if(glyco[j][10][0] && glyco[j][11][0]) ++supported;
  const size_t slots[6]={0,1,2,3,4,9};
  for(size_t k=0;k<6;++k) {
   size_t source=glyco[j][slots[k]][0],match=0,where=0;
   for(size_t i=0;i<coordinates.size();++i)
    if(input[source].coords[0]==coordinates[i][0] && input[source].coords[1]==coordinates[i][1] && input[source].coords[2]==coordinates[i][2]) {++match;where=i;}
   if(match!=1) throw std::runtime_error("non-bijective input atom correspondence");
   glyco[j][slots[k]][0]=where;
  }
 }
 conf initial=m.get_initial_conf();
 m.set(initial);
 double native=m.eval_chi(1.0,2.0),maximum=native;
 int changed=-1;
 conf chosen=initial;
 for(size_t i=0;i<initial.ligands[0].torsions.size();++i) {
  conf trial=initial;trial.ligands[0].torsions[i]+=1.57079632679489661923;
  m.set(trial);double value=m.eval_chi(1.0,2.0);
  if(value>maximum) {maximum=value;chosen=trial;changed=i;}
 }
 m.set(chosen);
 double zero=m.eval_chi(0.0,2.0),one=m.eval_chi(1.0,2.0),two=m.eval_chi(2.0,2.0);
 std::cout<<std::setprecision(17)<<"{\"glycosidic_links\":"<<glyco.size()<<",\"fully_supported_links\":"<<supported<<",\"initial_chi_kcal_mol\":"<<native<<",\"torsion_perturbed\":"<<changed<<",\"coefficient_0\":"<<zero<<",\"coefficient_1\":"<<one<<",\"coefficient_2\":"<<two<<"}\n";
 return 0;
}
