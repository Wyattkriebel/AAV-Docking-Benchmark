#include <iostream>
#include <iomanip>
#include <stdexcept>
#include "model.h"
#include "parse_pdbqt.h"
#include "everything.h"
#include "weighted_terms.h"
#include "naive_non_cache.h"
#include "precalculate.h"
int main(int argc,char** argv) {
 if(argc!=3) return 2;
 model m=parse_receptor_pdbqt(path(argv[1]));
 m.append(parse_ligand_pdbqt(path(argv[2])));
 vecv xyz=m.get_flexible_coords();
 std::vector<parsed_atom> input=liginfo_return();
 std::vector<std::vector<size_t*> > glyco=glycan_info_func();
 size_t supported=0;
 for(size_t j=0;j<glyco.size();++j) {
  if(glyco[j][10][0] && glyco[j][11][0]) ++supported;
  const size_t slots[6]={0,1,2,3,4,9};
  for(size_t k=0;k<6;++k) {
   size_t source=glyco[j][slots[k]][0],match=0,where=0;
   for(size_t i=0;i<xyz.size();++i)
    if(input[source].coords[0]==xyz[i][0] && input[source].coords[1]==xyz[i][1] && input[source].coords[2]==xyz[i][2]) {++match;where=i;}
   if(match!=1) throw std::runtime_error("non-bijective ligand atom correspondence");
   glyco[j][slots[k]][0]=where;
  }
 }
 everything terms;
 flv weights;
 weights.push_back(-0.035579);weights.push_back(-0.005156);weights.push_back(0.840245);
 weights.push_back(-0.035069);weights.push_back(-0.587439);weights.push_back(5*0.05846/0.1-1);
 weighted_terms scoring(&terms,weights);
 precalculate prec(scoring);
 naive_non_cache receptor(&prec);
 const vec authentic(1000,1000,1000);
 conf c=m.get_initial_conf();
 double intra=m.eval_intramolecular(prec,authentic,c);
 double base0=m.eval_adjusted(scoring,prec,receptor,authentic,c,0.0);
 double score_only=m.eval_adjusted(scoring,prec,receptor,authentic,c,intra);
 double chi0=m.eval_chi(0.0,2.0),chi1=m.eval_chi(1.0,2.0);
 double normalization=scoring.conf_independent(m,1.0)-scoring.conf_independent(m,0.0);
 std::cout<<std::setprecision(17)<<"{\"glycosidic_links\":"<<glyco.size()<<",\"fully_supported_links\":"<<supported<<",\"intramolecular_reference\":0,\"base_shared0\":"<<base0<<",\"chi_coefficient0\":"<<chi0<<",\"chi_coefficient1\":"<<chi1<<",\"combined_shared0\":"<<(base0+chi1)<<",\"original_score_only_base\":"<<score_only<<",\"pose_intramolecular_energy\":"<<intra<<",\"torsion_normalization_scale\":"<<normalization<<",\"nonbonded_cutoff_A\":"<<scoring.cutoff()<<"}\n";
 return 0;
}
