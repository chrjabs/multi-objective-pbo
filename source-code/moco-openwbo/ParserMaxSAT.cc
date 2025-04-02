#include "ParserMaxSAT.h"

using Glucose::Solver;
namespace openwbo {
  void newSATVariable(Solver *S) {
#ifdef SIMP
    ((NSPACE::SimpSolver *)S)->newVar();
#else
    S->newVar();
#endif
  }
  
  int get_new_Lit_id(Solver * S){
    Lit p = mkLit(S->nVars(), false);
    newSATVariable(S);
    return var(p);
  }


  void parseMOCNF(char * fname, Solver *S, tmp_kp_t * tmp, std::map<std::string, int> name2id, std::map<int, bool> id2sign, Lit relax_var)
  {
    gzFile input_stream = gzopen(fname, "rb");
    if (input_stream == NULL)
      printf("c ERROR! Could not open file: %s\n", fname),
	printf("s UNKNOWN\n"), exit(_ERROR_);
    
    StreamBuffer in(input_stream);
    //         
      vec<Lit> lits;
  
      //   int nObj = -1;
      std::map<int, int> fid2hid; //from id in file to id in maxsat_formula (here)
      std::map<int, std::string> fid2x; //id assigned to variables x in the file

      std::vector<int> tmp_vals; //tmp
      int nls = 0; //count how many were read
    
      tmp->haszs = false;
    
      //     printf("c Solver before parseMOCNF\n");
      //     S->my_print();
      
      for (;;) {
	skipWhitespace(in);
	if (*in == EOF)
	  break;
	else if (*in == 'p') {
	  if (eagerMatch(in, "p mocnf")) {
	    parseInt(in); // Variables
	    parseInt(in); // Clauses
	    parseInt(in); // Objectives
	    parseInt(in); // input vars
	    parseInt(in); // output vars
        
	  }else
	    printf("c PARSE ERROR! Unexpected char: %c\n", *in),
	      printf("s UNKNOWN\n"), exit(_ERROR_);
	} else if (*in == 'c' || *in == 'p')
	  skipLine(in);
	else if (*in == 'b'){
	  ++in;
	  parseInt(in); // obj id
         
	  readLineVals(in, tmp_vals);
         
	  tmp->nbase = tmp_vals.size();
	  tmp->ncoeffs = tmp->nbase + 1;
	  //          printf("base size: %d\n", tmp->nbase);
	  //          printf("ncoeffs size: %d\n", tmp->ncoeffs);
         
	  tmp->base = (uint64_t *) malloc(tmp->nbase * sizeof(uint64_t));
	  tmp->mxb_coeffs = (uint64_t *) malloc(tmp->ncoeffs * sizeof(uint64_t));
         
	  tmp->ls = (int **) malloc(tmp->ncoeffs * sizeof(int *));
	  tmp->ls_szs = (uint64_t *) malloc(tmp->ncoeffs * sizeof(uint64_t));
         
	  for(int i = 0; i < tmp->ncoeffs; i++){
	    tmp->ls[i] = NULL; 
	  }

	  tmp->mxb_coeffs[0] = 1;
	  tmp->ls_szs[0] = 0;
	  for(int i = 0; i < tmp->nbase; i++){
	    tmp->base[i] = tmp_vals[i];
	    tmp->mxb_coeffs[i+1] = tmp->mxb_coeffs[i] * tmp->base[i];
	    tmp->ls_szs[i+1] = 0;
	  }
         
	  tmp->nvars = 0;
	  tmp->nclauses = 0;
         
	} else if (*in == 'l' || *in == 'z'){
	  //          printf("l or z\n");
	  if(*in == 'z') tmp->haszs = true;
         
	  ++in;
	  parseInt(in); // obj id
	  parseInt(in); // coeficiente
	  //          printf("di, coeff: %d %d\n", di, coeff);
         
	  readLineVals(in, tmp_vals);
	  //          printf("coeff: %lu\n", coeff);
	  //          printf("tmp->mxb_coeffs[nls]: %lu\n", tmp->mxb_coeffs[nls]);
	  //          printf("nls: %d, tmp->base[nls]: %lu\n", nls, tmp->base[nls]);
         
	  //          assert(coeff == tmp->mxb_coeffs[nls]);
	  int size = tmp_vals.size();
	  //          printf("ls_ids: %d\n", size);
	  tmp->ls_szs[nls] = size;
	  tmp->ls[nls] = (int *) malloc(size * sizeof(int));
	  for(int j = 0; j < size; j++){
	    tmp->ls[nls][j] = tmp_vals[j];
	    //                 printf(" %d", tmp->ls[nls][j]);
	  }
	  //         printf("\n");
        
	  nls++;


	}else if (*in == 'x'){
	  //          printf("x\n");
	  ++in; ++in;
	  static char word[10];
	  int sz, id;
	  parseOriginalVar(in, word, &sz);
	  id = parseInt(in);
	  fid2x[id] = std::string(word);
         
	  //          printf("%s: %d\n", word, id);
	}
	else {
        
	  if(fid2hid.size() == 0){
	    //print
	    /*
	      printf("oname solverid (size: %d)\n", name2id.size());
	      for(std::map<std::string, int>::iterator it = name2id.begin(); it != name2id.end(); it++){
	      printf("%s -> %d\n", it->first.c_str(), it->second);
	      }
	    */

	int oid, t;
	  for(std::map<int, std::string>::iterator it = fid2x.begin(); it != fid2x.end(); it++){
	    fid2hid[it->first] = name2id.at(it->second);
	    //TODO: Ter um dicionario de sinais, e aqui verificar qual o sinal destes
	    //lits (da funcao objectivo) e usá-los na construcao das clausulas...
	    oid = fid2hid.at(it->first);
	    //se a var original tem coef negativo, criar uma var aux para inverter o sinal
	    if(!id2sign.at(oid)){
	      t = get_new_Lit_id(S);
	      //                     printf("invert x: ~%d -> %d\n", oid, t);
	      lits.clear();
	      //t => ~oid
	      if(relax_var != lit_Undef)
		lits.push(relax_var);
	      lits.push(~mkLit(t)); lits.push(~mkLit(oid)); S->addClause(lits);
	      lits.clear();
	      //~oid => t
	      if(relax_var != lit_Undef)
		lits.push(relax_var);
	      lits.push(mkLit(t)); lits.push(mkLit(oid)); S->addClause(lits);
	      fid2hid[it->first] = t;
                    
	      tmp->nclauses+=2;
	    }
	  }
	  //get new ids for ls

	  int l, k;
	  int size;
	  for(int i = 0; i < tmp->ncoeffs; i++){
	    size = tmp->ls_szs[i];
	    for(int j = 0; j < size; j++){
	      l = tmp->ls[i][j];
	      if(fid2hid.find(l) == fid2hid.end()){ //if key does not exist in fid2hid
		//                         printf("%d is not in fid2hid\n", l);
		k = get_new_Lit_id(S);
		fid2hid[l] = k;
	      }else{
		k = fid2hid.at(l);
	      }
	      tmp->ls[i][j] = k; // fix ls -- use solver id instead of file id
	    }
	  }
	  //get new ids for zs
	  //TODO
        
        }

        //create clauses (get new ids for others)
        //TODO
        readLineVals(in, tmp_vals);
        lits.clear();
        int k;
        Lit l;
	//         printf("----------- KP clauses --------------------\n");
        for(uint i = 0; i < tmp_vals.size(); i++){
	  int v = abs(tmp_vals[i]);
	  if(fid2hid.find(v) == fid2hid.end()){
	    //                 printf("%d is not in fid2hid\n", v);
	    k = get_new_Lit_id(S);
	    //                 printf("new lit id: %d\n", k);
	    fid2hid[v] = k;
	  }
	  k = fid2hid[v];
	  l = (tmp_vals[i] > 0) ? mkLit(k) : ~mkLit(k);
	  lits.push(l);
	  //             printf(" %c%d[%d]", (tmp_vals[i] > 0) ? ' ':'-', v, k);
        }
	//         printf("\n");
        
        S->addClause(lits);
        tmp->nclauses++;
        tmp->nvars = fid2hid.size();
        
	//         S->my_print();

      
      }
    }
  
    // ------------ prints -----------------------------------------------------------------
    /*
      printf("\nc base (%d):\n", tmp->nbase);
      for(int i = 0; i < tmp->nbase; i++) printf(" %d", tmp->base[i]);
    
      char s;
      printf("\n\nc coeffs n ls (%d):\n", tmp->ncoeffs);
      for(int i = 0; i < tmp->ncoeffs; i++){
      printf("c [%d]", tmp->mxb_coeffs[i]);
      for(int j = 0; j < tmp->ls_szs[i]; j++){
      //             if(id2sign.find(tmp->ls[i][j]) != id2sign.end())
      //                 s = (id2sign.at(tmp->ls[i][j])) ? '+' : '-';
      //             printf(" %c%d", s, tmp->ls[i][j]);
      printf(" %d", tmp->ls[i][j]);
      }
      printf("\n");
      }*/
    /*
      printf("\nfileid solverid\n");
      for(std::map<int, int>::iterator it = fid2hid.begin(); it != fid2hid.end(); it++){
printf("%d -> %d\n", it->first, it->second);
    }
    */
//     for(int i = 0; i < nls; i++) free(ls_ids[i]);
//     free(ls_ids);
//     free(ls_szs);

//     S->my_print();
    
    // -------------------------------------------------------------------------------------
    
    gzclose(input_stream);
}
}
