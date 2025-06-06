#include <vector>
#include <queue>
#include <algorithm>
#include <cassert>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <set>
#include <sstream>
#include <stack>


#include "preprocessor.hpp"
#include "trace.hpp"
#include "global.hpp"
#include "utility.hpp"
#include "timer.hpp"
#include "log.hpp"
#include "AMSLEX.hpp"
#include "satlikeinterface.hpp"

#define F first
#define S second

using namespace std;
namespace maxPreprocessor {

void Preprocessor::init() {
	originalVars = pi.vars;
	originalClauses = pi.clauses.size();
	BIGIt = 1;
	BVElocalGrow = 0;
	BVEglobalGrow = 0;

	redSatSolverCalls = 0;

	bestCost = HARDWEIGHT;

	doneUnhiding = false;
	randGen.seed(123);
	flePos = 0;
	fleActiveTechniques = 0;
	rLog.weightRange = pi.getWeightSums();
	rLog.initialWeightRange = rLog.weightRange;

	for (int i=0; i<pi.objectives; ++i) {
		stats["original_weight_sum"+std::to_string(i)]=rLog.weightRange[i];
	}
	stats["original_clauses"]=pi.clauses.size();
	stats["original_variables"]=pi.vars;
}

Preprocessor::Preprocessor(const vector<vector<int> >& clauses, const vector<vector<uint64_t> >& weights, uint64_t topWeight)
	: pi(clauses, weights, topWeight), satSolver(0), amsLex(pi) {
		init();
}
Preprocessor::Preprocessor(const vector<vector<int> >& clauses, const vector<uint64_t>& weights, uint64_t topWeight)
	: pi(clauses, weights, topWeight), satSolver(0), amsLex(pi) {
		init();
}
Preprocessor::Preprocessor(const vector<vector<int> >& clauses, const vector<pair<uint64_t, uint64_t> >& weights, uint64_t topWeight)
	: pi(clauses, weights, topWeight), satSolver(0), amsLex(pi) {
		init();
}


void Preprocessor::prepareSatSolver() {
	// For now use always glucose 3... and create it every time...
	if (satSolver) {
		delete satSolver;
	}
	satSolver = (SATSolver*) new Glucose3();

	for (unsigned c = 0; c<pi.clauses.size(); ++c) {
		if (pi.isClauseRemoved(c)) continue;
		if (!pi.clauses[c].isHard()) continue;

		satSolver->addClause(pi.clauses[c].lit);
	}
}

bool Preprocessor::isTautology(const Clause& clause) const {
	for (unsigned i = 1; i < clause.lit.size(); i++) {
		if (litNegation(clause.lit[i]) == clause.lit[i - 1]) {
			return true;
		}
	}
	return false;
}

// Returns number of clauses removed
int Preprocessor::setVariable(int var, bool value) {
	int removed = 0;
	trace.setVar(var, value);
	vector<int>& satClauses = (value == true) ? pi.litClauses[posLit(var)] : pi.litClauses[negLit(var)];
	vector<int>& notSatClauses = (value == true) ? pi.litClauses[negLit(var)] : pi.litClauses[posLit(var)];
	for (int c : satClauses) {
		pi.removeClause(c);
		removed++;
	}
	for (int c : notSatClauses) {
		if (value == true) pi.removeLiteralFromClause(negLit(var), c);
		else pi.removeLiteralFromClause(posLit(var), c);
	}
	return removed;
}

// This is called only in the beginning since no tautologies are added
void Preprocessor::removeTautologies() {
	int found = 0;

	for (unsigned i = 0; i < pi.clauses.size(); i++) {
		if (isTautology(pi.clauses[i])) {
			found++;
			pi.removeClause(i);
		}
	}

	log(found, " tautologies removed");
}

int Preprocessor::eliminateReduntantLabels() {
	int fnd = 0;

	for (int objective=0; objective<pi.objectives; ++objective) {
		map<int64_t, map<int, vector<int> > > labels;
		vector<int> ls;
		for (int var = 0; var < pi.vars; var++) {
			if (pi.labelIndexMask(var) == (1<<objective)) {
				if (pi.labelPolarity(var, objective) == VAR_TRUE) {
					if (pi.litClauses[posLit(var)].size() == 1 && pi.litClauses[negLit(var)].size() == 1) {
						int c = pi.litClauses[negLit(var)][0];
						for (int lit : pi.clauses[c].lit) {
							labels[pi.labelWeight(var, objective)][lit].push_back(posLit(var));
						}
						ls.push_back(posLit(var));
					}
				}
				else {
					if (pi.litClauses[negLit(var)].size() == 1 && pi.litClauses[posLit(var)].size() == 1) {
						int c = pi.litClauses[posLit(var)][0];
						for (int lit : pi.clauses[c].lit) {
							labels[pi.labelWeight(var, objective)][lit].push_back(negLit(var));
						}
						ls.push_back(negLit(var));
					}
				}
			}
		}
		vector<int> matched(pi.vars);
		//priorize matching ternary clauses
		for (int lb : ls) {
			if (matched[litVariable(lb)]) continue;
			int c = pi.litClauses[litNegation(lb)][0];
			if (pi.clauses[c].lit.size() != 3) continue;
			vector<int> lits;
			for (int lit : pi.clauses[c].lit) {
				if (!pi.isLabelVar(litVariable(lit))) lits.push_back(litNegation(lit));
			}
			if (lits.size() != 2) continue;
			int tl = lits[0];
			if (pi.litClauses[lits[1]].size() < pi.litClauses[tl].size()) tl = lits[1];
			for (int c2 : pi.litClauses[tl]) {
				if (pi.clauses[c2].lit.size() != 3) continue;
				int lb2 = -1;
				vector<int> lits2;
				for (int lit : pi.clauses[c2].lit) {
					if (!pi.isLabelVar(litVariable(lit))) lits2.push_back(lit);
					else if (pi.labelIndexMask(lit) != (1<<objective)) break;
					else lb2 = litNegation(lit);
				}
				if (lits2.size() != 2) continue;
				if (lb2 == -1) continue;
				if (lb2 == lb) continue;
				if (pi.labelWeight(litVariable(lb), objective) != pi.labelWeight(litVariable(lb2), objective)) continue;
				if (matched[litVariable(lb2)]) continue;
				if (pi.litClauses[litNegation(lb2)].size() > 1) continue;
				if ((lits[0] == lits2[0] && lits[1] == lits2[1]) || (lits[1] == lits2[0] && lits[0] == lits2[1])) {
					fnd++;
					matched[litVariable(lb)] = 1;
					matched[litVariable(lb2)] = 1;
					pi.removeLiteralFromClause(litNegation(lb2), c2);
					pi.addLiteralToClause(litNegation(lb), c2);
					pi.removeClause(pi.litClauses[lb2][0]);
					assert(pi.isVarRemoved(litVariable(lb2)));
					trace.labelEliminate(lb, lb2, lits2[0]);
					trace.setVar(litVariable(lb2), litValue(lb2));
					pi.unLabel(litVariable(lb2), objective);
					break;
				}
			}
		}
		//match greedily
		for (unsigned i = 0; i < ls.size(); i++) {
			if (matched[litVariable(ls[i])]) continue;
			int c1 = pi.litClauses[litNegation(ls[i])][0];
			for (int lit : pi.clauses[c1].lit) {
				if (lit == litNegation(ls[i])) continue;
				for (int l2 : labels[pi.labelWeight(litVariable(ls[i]), objective)][litNegation(lit)]) {
					if (matched[litVariable(l2)]) continue;
					int c2 = pi.litClauses[litNegation(l2)][0];
					bool taut = false;
					int tautli = 0;
					unsigned j2 = 0;
					for (unsigned j = 0; j < pi.clauses[c1].lit.size(); j++) {
						while (j2 < pi.clauses[c2].lit.size() && pi.clauses[c2].lit[j2] < pi.clauses[c1].lit[j]) {
							if (pi.clauses[c2].lit[j2] == litNegation(pi.clauses[c1].lit[j])) {
								tautli = pi.clauses[c2].lit[j2];
								taut = true;
							}
							j2++;
						}
						if (j2 < pi.clauses[c2].lit.size() && pi.clauses[c2].lit[j2] == litNegation(pi.clauses[c1].lit[j])) {
							tautli = pi.clauses[c2].lit[j2];
							taut = true;
						}
						if (taut) break;
					}
					assert(taut == true);
					fnd++;
					matched[litVariable(ls[i])] = 1;
					matched[litVariable(l2)] = 1;
					pi.removeLiteralFromClause(litNegation(l2), c2);
					pi.addLiteralToClause(litNegation(ls[i]), c2);
					pi.removeClause(pi.litClauses[l2][0]);
					assert(pi.isVarRemoved(litVariable(l2)));
					trace.labelEliminate(ls[i], l2, tautli);
					trace.setVar(litVariable(l2), litValue(l2));
					pi.unLabel(litVariable(l2), objective);
					break;
				}
				if (matched[litVariable(ls[i])]) break;
			}
		}
	}
	return fnd;
}

// This is called only in the beginning
void Preprocessor::identifyLabels() {
	int found = 0;
	int eliminated = 0;

	// Find unit soft clauses where the negation of the literal occurs only in hard clauses
	for (int lit = 0; lit < pi.vars*2; lit++) {
		// if (pi.isLabelVar(litVariable(lit))) continue;
// 		if (pi.litClauses[lit].size()>1) continue;
		bool f = false;
		vector<int> toRemove;
		for (unsigned i=0; i<pi.litClauses[lit].size(); ++i) {
			int c=pi.litClauses[lit][i];
			if (!pi.clauses[c].isHard() && pi.clauses[c].lit.size() == 1) {
				if (!f) {
					f = true;
					swap(pi.litClauses[lit][i], pi.litClauses[lit][0]);
				} else {
					pi.pourAllWeight(c, pi.litClauses[lit][0]);
					toRemove.push_back(c);
				}
			}
		}
		if (!f) continue;
		// check unit negations...
		f = false;
		for (int c : pi.litClauses[litNegation(lit)]) {
			if (!pi.clauses[c].isHard() && pi.clauses[c].lit.size() == 1) {
				trace.removeWeight(pi.substractWeights(c, pi.litClauses[lit][0]));


				if (pi.clauses[c].isWeightless()) {
					toRemove.push_back(c);
				}

				if (pi.clauses[pi.litClauses[lit][0]].isWeightless()) {
					toRemove.push_back(pi.litClauses[lit][0]);
					f = true;
					break;
				}
			}
		}
		eliminated += toRemove.size();
		for (int c : toRemove) {
			pi.removeClause(c);
		}
		if (f) continue;

		if (litValue(lit) == true) {
			for (int objective=0; objective<pi.objectives; ++objective) {
				if (pi.clauses[pi.litClauses[lit][0]].weight(objective))
					pi.mkLabel(litVariable(lit), objective, VAR_TRUE);
			}
		} else {
			for (int objective=0; objective<pi.objectives; ++objective) {
				if (pi.clauses[pi.litClauses[lit][0]].weight(objective))
					pi.mkLabel(litVariable(lit), objective, VAR_FALSE);
			}
		}
		found++;
	}

	log(found, " labels identified");
	log(eliminated, " soft unit clauses eliminated");
}

// This is called only in the beginning
void Preprocessor::createLabels() {
	int added = 0;

	// Create labels for every soft clause that does not have a label yet
	for (unsigned i = 0; i < pi.clauses.size(); i++) {
		if (!pi.clauses[i].isHard() && !pi.isClauseRemoved(i) && !pi.isLabelClause(i)) {
			int nv = pi.addVar();
			pi.addLiteralToClause(posLit(nv), i);
			pi.addClause({negLit(nv)}, pi.clauses[i].weights);
			for (int objective=0; objective<pi.objectives; ++objective) {
				if (pi.clauses[i].weight(objective))
					pi.mkLabel(nv, objective, VAR_FALSE);
			}
			pi.clauses[i].makeHard();
			added++;
		}
	}

	log(added, " labels added");
}

int Preprocessor::removeEmptyClauses() {
	int removed = 0;
	vector<int> src;
	for (unsigned i = 0; i < pi.clauses.size(); i++) {
		if (!pi.isClauseRemoved(i)) {
			if (pi.clauses[i].lit.size() == 0) {
				if (pi.clauses[i].isHard()) {

				}
				else {
					src.push_back(i);
				}
			}
		}
	}
	for (int c : src) {
		trace.removeWeight(pi.clauses[c].weights);
		pi.removeClause(c);
		removed++;
	}
	log(removed, " empty clauses removed");
	return removed;
}

int Preprocessor::tryUP(int lit) {
	for (int c : pi.litClauses[lit]) {
		if (pi.clauses[c].lit.size() == 1 && pi.clauses[c].isHard()) {
			if (pi.isLabelVar(litVariable(lit))) {
				rLog.removeLabel(1);
			}	else {
				rLog.removeVariable(1);
			}
			int rmClauses = setVariable(litVariable(lit), litValue(lit));
			rLog.removeClause(rmClauses);
			return rmClauses;
		}
	}
	return 0;
}

int Preprocessor::tryUPAll() {
	int removed = 0;
	for (unsigned c = 0; c < pi.clauses.size(); c++) {
		if (pi.clauses[c].lit.size() == 1 && pi.clauses[c].isHard() && !pi.isClauseRemoved(c)) {
			int lit = pi.clauses[c].lit[0];
			if (pi.isLabelVar(litVariable(lit))) {
				rLog.removeLabel(1);
			}	else {
				rLog.removeVariable(1);
			}
			int rmClauses = setVariable(litVariable(lit), litValue(lit));
			rLog.removeClause(rmClauses);
			removed += rmClauses;
			if (!rLog.requestTime(Log::Technique::UP)) break;
		}
	}
	return removed;
}

int Preprocessor::doUP() {
	rLog.startTechnique(Log::Technique::UP);
	int removed = 0;
	if (!rLog.requestTime(Log::Technique::UP)) {
		rLog.stopTechnique(Log::Technique::UP);
		return 0;
	}
	vector<int> checkLit = pi.tl.getTouchedLiterals("UP");
	if ((int)checkLit.size() > pi.vars) {
		removed = tryUPAll();
	}
	else {
		for (int lit : checkLit) {
			if (!rLog.requestTime(Log::Technique::UP)) break;
			removed += tryUP(lit);
		}
	}

	log(removed, " clauses removed by UP");
	rLog.stopTechnique(Log::Technique::UP);
	return removed;
}

void Preprocessor::doUP2() {
	rLog.startTechnique(Log::Technique::UP);
	for (int lit = 0; lit < 2*pi.vars; lit++) {
		assert(tryUP(lit) == 0);
	}
	rLog.stopTechnique(Log::Technique::UP);
}

int Preprocessor::removeDuplicateClauses() {
	int removed = 0;
	vector<int> checkLit = pi.tl.getModLiterals("DPCLRM");
	vector<pair<uint64_t, int> > has;
	if ((int)checkLit.size() > pi.vars/2) { // magic constant
		for (int c = 0; c < (int)pi.clauses.size(); c++) {
			if (!pi.isClauseRemoved(c)) {
				if (pi.clauses[c].lit.size() == 0) {
					if (!pi.clauses[c].isHard()) {
						trace.removeWeight(pi.clauses[c].weights);
						pi.removeClause(c);
					}
				}
				else {
					uint64_t h = 0;
					for (int l : pi.clauses[c].lit) {
						h *= polyHashMul;
						h += (uint64_t)(l + 1);
					}
					has.push_back({h, c});
				}
			}
		}
		sort(has.begin(), has.end());
		bool hard = false;
		vector<int> toRemove;
		for (unsigned i = 1; i < has.size(); i++) {
			if (has[i].F != has[i - 1].F) {
				hard = false;
				continue;
			}
			if (pi.clauses[has[i].S].lit != pi.clauses[has[i - 1].S].lit) {
				hard = false;
				continue;
			}

			if (hard || pi.clauses[has[i-1].S].isHard()) {
				hard = true;
				toRemove.push_back(has[i].S);
			} else if (pi.clauses[has[i].S].isHard()) {
				hard = true;
				toRemove.push_back(has[i-1].S);
			} else {
				pi.pourAllWeight(has[i-1].S, has[i].S);
				toRemove.push_back(has[i-1].S);
			}
		}
		for (int c : toRemove) assert(!pi.isClauseRemoved(c));
		for (int c : toRemove) {
			if (pi.isLabelClause(c)) { // since after labelling there are no duplicate softs, this implies we have a unit hard clause and we can harden it
				removed += setVariable(litVariable(pi.clauses[c].lit[0]), litValue(pi.clauses[c].lit[0]));
			} else if (!pi.isClauseRemoved(c)) { // hardening can remove some clauses marked as to remove...
				removed++;
				pi.removeClause(c);
			}
		}
	}
	else {
		for (int lit : checkLit) {
			has.clear();
			for (int c : pi.litClauses[lit]) {
				uint64_t h = 0;
				for (int l : pi.clauses[c].lit) {
					h *= polyHashMul;
					h += (uint64_t)(l + 1);
				}
				has.push_back({h, c});
			}
			sort(has.begin(), has.end());
			bool hard = false;
			vector<int> toRemove;
			for (unsigned i = 1; i < has.size(); i++) {
				if (has[i].F != has[i - 1].F) {
					hard = false;
					continue;
				}
				if (pi.clauses[has[i].S].lit != pi.clauses[has[i - 1].S].lit) {
					hard = false;
					continue;
				}

				if (hard || pi.clauses[has[i-1].S].isHard()) {
					hard = true;
					toRemove.push_back(has[i].S);
				} else if (pi.clauses[has[i].S].isHard()) {
					hard = true;
					toRemove.push_back(has[i-1].S);
				} else {
					pi.pourAllWeight(has[i-1].S, has[i].S);
					toRemove.push_back(has[i-1].S);
				}
			}
			for (int c : toRemove) assert(!pi.isClauseRemoved(c));
			for (int c : toRemove) {
				if (pi.isLabelClause(c)) { // since after labelling there is not duplicate softs, this implies we have unit hard clause and we can harden it
					removed+=setVariable(litVariable(pi.clauses[c].lit[0]), litValue(pi.clauses[c].lit[0]));
				} else if (!pi.isClauseRemoved(c)) { // hardening can remove some clauses marked as to remove...
					removed++;
					pi.removeClause(c);
				}
			}
		}
	}
	return removed;
}

#include "modelsearch.cpp"

#include "SE.cpp"
#include "BVE.cpp"
#include "SSR.cpp"
#include "BCE.cpp"
#include "SLE.cpp"
#include "BCR.cpp"
#include "SIE.cpp"
#include "BIG.cpp"
#include "BVA.cpp"
#include "GSLE.cpp"
#include "FLP.cpp"
#include "LS.cpp"
#include "AM1.cpp"
#include "TMS.cpp"
#include "RED.cpp"
#include "HARD.cpp"
#include "FLE.cpp"

PreprocessedInstance Preprocessor::getPreprocessedInstance(bool addRemovedWeight, bool sortLabelsFrequency) {
	PreprocessedInstance ret;
	for (unsigned i = 0; i < pi.clauses.size(); i++) {
		if (!pi.isClauseRemoved(i) && pi.clauses[i].isHard()) {
			ret.clauses.push_back(pi.clauses[i].lit);
			if (pi.objectives>1) ret.weightsv.push_back({});
			else                 ret.weights.push_back(HARDWEIGHT);
		}
	}
	for (int var = 0; var < pi.vars; var++) {
		if (pi.litClauses[posLit(var)].size() > 0 || pi.litClauses[negLit(var)].size() > 0) {
			if (pi.litClauses[negLit(var)].size() == 0) {
				trace.setVar(var, true);
			} else if (pi.litClauses[posLit(var)].size() == 0) {
				trace.setVar(var, false);
			} else {
				if (pi.isLitLabel(posLit(var))) {
					assert(pi.clauses[pi.litClauses[posLit(var)][0]].lit.size() == 1);
					assert(!pi.clauses[pi.litClauses[posLit(var)][0]].isHard());
					ret.labels.push_back({posLit(var), pi.clauses[pi.litClauses[posLit(var)][0]].weights});
				}
				if (pi.isLitLabel(negLit(var))) {
					assert(pi.clauses[pi.litClauses[negLit(var)][0]].lit.size() == 1);
					assert(!pi.clauses[pi.litClauses[negLit(var)][0]].isHard());
					ret.labels.push_back({negLit(var), pi.clauses[pi.litClauses[negLit(var)][0]].weights});
				}
			}
		}
	}

	if (addRemovedWeight && trace.removedWeight.size()) {
	 	int nVar = pi.getExcessVar();
	 	ret.labels.push_back({negLit(nVar), trace.removedWeight});
	 	ret.clauses.push_back({posLit(nVar)});
	 	if (pi.objectives>1) ret.weightsv.push_back({});
	 	else                 ret.weights.push_back(HARDWEIGHT);
	}

	if (sortLabelsFrequency) {
		auto cmp = [&](const pair<int, vector<uint64_t> >& a, const pair<int, vector<uint64_t> >&	 b) {
			if (litVariable(a.F) >= pi.vars || litVariable(b.F) >= pi.vars) return a.F > b.F;
			return pi.litClauses[a.F].size() + pi.litClauses[litNegation(a.F)].size() > pi.litClauses[b.F].size() + pi.litClauses[litNegation(b.F)].size();
		};
		sort(ret.labels.begin(), ret.labels.end(), cmp);
	}
	for (int i = 0; i < (int)ret.labels.size(); i++) {
		ret.clauses.push_back({ret.labels[i].F});
		if (pi.objectives>1) ret.weightsv.push_back(ret.labels[i].S);
		else                 ret.weights.push_back(ret.labels[i].S[0]);
	}
	return ret;
}

bool Preprocessor::validTechniques(string techniques) const {
	int sb = 0;
	string vt = "buvsrilceagphtmGSQTVdDMLHURP";
	for (int i = 0; i < (int)techniques.size(); i++) {
		if(techniques[i] == '[') {
			sb++;
		}
		else if(techniques[i] == ']') {
			sb--;
			if (sb < 0) return false;
			if (techniques[i - 1] == '[') return false;
		}
		else {
			bool f = false;
			for (char c : vt) {
				if (techniques[i] == c) {
					f = true;
					break;
				}
			}
			if (!f) return false;
		}
	}
	if (sb != 0) return false;
	return true;
}

bool Preprocessor::validPreTechniques(string techniques) const {
	int sb = 0;
	string vt = "busU";
	for (int i = 0; i < (int)techniques.size(); i++) {
		if(techniques[i] == '[') {
			sb++;
		}
		else if(techniques[i] == ']') {
			sb--;
			if (sb < 0) return false;
			if (techniques[i - 1] == '[') return false;
		}
		else {
			bool f = false;
			for (char c : vt) {
				if (techniques[i] == c) {
					f = true;
					break;
				}
			}
			if (!f) return false;
		}
	}
	if (sb != 0) return false;
	return true;
}

//b = BCE, u = UP, v = BVE, s = SE, r = SSR, l = SLE, c = BCR, i = SIE, e = EE, a = BVA, g = GSLE, p = FLP, h = UH, t = structure labeling, m = AM1, T = TMS, d=LRED, D=CRED, H=HARD

int Preprocessor::doPreprocess(const string& techniques, int l, int r, bool debug, bool topLevel) {
	int found = 0;
	if (l == r) {
		char tc = techniques[l];
		if (tc == 'b') {
			while (int f = doBCE()) {
				found += f;
			}
			if (debug) doBCE2();
		}
		else if (tc == 'u') {
			while (int f = doUP()) {
				found += f;
			}
			if (debug) doUP2();
		}
		else if (tc == 'v') {
			while (int f = doBVE()) {
				found += f;
			}
			if (debug) doBVE2();
		}
		else if (tc == 's') {
			while (int f = doSE()) {
				found += f;
			}
			if (debug) doSE2();
		}
		else if (tc == 'r') {
			while (int f = doSSR()) {
				found += f;
			}
			if (debug) doSSR2();
		}
		else if (tc == 'l') {
			while (int f = doSLE()) {
				found += f;
			}
			if (debug) doSLE2();
		}
		else if (tc == 'c') {
			while (int f = doBCR()) {
				found += f;
			}
			if (debug) doBCR2();
		}
		else if (tc == 'i') {
			while (int f = doSIE()) {
				found += f;
			}
			if (debug) doSIE2();
		}
		else if (tc == 'e') {
			while (int f = doBIG(false)) {
				found += f;
			}
			if (debug) doBIG2(false);
		}
		else if (tc == 'a') {
			while (int f = doBVA()) {
				found += f;
			}
			if (debug) doBVA2();
		}
		else if (tc == 'g') {
			while (int f = doGSLE()) {
				found += f;
			}
			if (debug) doGSLE2();
		}
		else if (tc == 'p') {
			while (int f = doFLP()) {
				found += f;
			}
			// FLP already checks every label
		}
		else if (tc == 'h') {
			while (int f = doBIG(true)) {
				found += f;
			}
			if (debug) doBIG2(true);
		}
		else if (tc == 't') {
			found += doLS();
		}
		else if (tc == 'm') {
			found += doAM1(1, 0, 0);
		}
		else if (tc == 'G') {
			found += doAM1(1, 0, 1);
		}
		else if (tc == 'S') {
			found += doAM1(1, 1, 1);
		}
		else if (tc == 'Q') {
			while (int f = doAM1(1, 1, 1)) {
				found += f;
			}
		}
		else if (tc == 'T') {
			found += doTMS();
		}
		else if (tc == 'V') {
			while (int f=doBBTMS(opt.BBTMS_maxVars)) {
				found += f;
			}
		}
		else if (tc == 'd') {
			found += doLabelRED();
		}
		else if (tc == 'D') {
			found += doClauseRED();
		}
		else if (tc == 'M') {
			found += doModelCuttingRED();
		}
		else if (tc == 'L') {
			found += doUPLitRED();
		}
		else if (tc == 'H') {
			found += doHARD();
		}
		else if (tc == 'U') {
			found += doFLE(0,1,0,1,0,0);
		}
		else if (tc == 'R') {
			found += doFLE(1,1,1,1,1,opt.FLE_redTechniques);
		}
		else {
			abort();
		}
		pi.tl.shrink(logLevel > 0);
	}
	else {
		int lp = l;
		int br = 0;
		for (int i = l; i <= r; i++) {
			if (techniques[i] == '[') {
				br++;
			}
			else if(techniques[i] == ']') {
				br--;
				if (br == 0) {
					while (int f = doPreprocess(techniques, lp + 1, i - 1, debug, false)) {
						found += f;
					}
					if (topLevel) {
						for (int j = lp + 1; j <= i - 1; j++) {
							if (techniques[j] != ']' && techniques[j] != '[' && (techniques[j] < '0' || techniques[j] > '9')) {
								bool neverAgain = true;
								for (int jj = i + 1; jj <= r; jj++) {
									if (techniques[jj] == techniques[j]) {
										neverAgain = false;
									}
								}
								if (neverAgain) {
									rLog.neverAgain(techniques[j]);
								}
							}
						}
					}
					lp = i + 1;
				}
			}
			else {
				if (br == 0) {
					assert(lp == i);
					found += doPreprocess(techniques, i, i, debug, false);
					if (topLevel) {
						bool neverAgain = true;
						for (int j = i + 1; j <= r; j++) {
							if (techniques[j] == techniques[i]) {
								neverAgain = false;
							}
						}
						if (neverAgain) {
							rLog.neverAgain(techniques[i]);
						}
					}
					lp = i + 1;
				}
			}
		}
		assert(br == 0);
	}
	return found;
}

void Preprocessor::preprocess(string techniques, double timeLimit, bool debug, bool BVEgate_, bool initialCall, bool matchLabels) {
	opt.BVEgate = BVEgate_;
	Timer preTime;
	preTime.start();
	log(originalVars, " variables, ", originalClauses, " clauses");

	print("c techniques ", techniques);
	log("techniques ", techniques);

	string preTechniques;
	for (unsigned i = 0; i < techniques.size(); i++) {
		if (techniques[i] == '#') {
			preTechniques = techniques.substr(0, i);
			techniques = techniques.substr(i + 1);
			break;
		}
	}

	if (!validTechniques(techniques) || !validPreTechniques(preTechniques)) {
		log("Invalid techniques");
		print("c Invalid techniques");
		abort();
	}

	if (initialCall) {
		removeTautologies();
		pi.tl.init(pi.vars);
		removeEmptyClauses();
	}

	int dpRm = removeDuplicateClauses();
	log(dpRm, " duplicate clauses removed");

	if (initialCall) {
		if (preTechniques.size() > 0) {
			preTime.stop(); // be sure that these 3 lines will be ran exactly once
			rLog.timePlan(timeLimit - preTime.getTime().count()*2, techniques);
			rLog.startTimer();
			doPreprocess(preTechniques, 0, (int)preTechniques.size() - 1, debug, false);
		}

		removeEmptyClauses();
		identifyLabels();
		createLabels();

		if (matchLabels) {
			int labelsMatched = eliminateReduntantLabels();
			rLog.labelsMatched += labelsMatched;
			log(labelsMatched, " labels matched");
		}

		dpRm = removeDuplicateClauses();
		log(dpRm, " duplicate clauses removed");

		pi.tl.init(pi.vars);
	}

	if (!initialCall || preTechniques.size() == 0) {
		preTime.stop();// here
		rLog.timePlan(timeLimit - preTime.getTime().count(), techniques);
		rLog.startTimer();
	}

	doPreprocess(techniques, 0, (int)techniques.size() - 1, debug, true);

	rLog.stopTimer();

	removeEmptyClauses();

	dpRm = removeDuplicateClauses();
	log(dpRm, " duplicate clauses removed");

	rLog.weightRange = pi.getWeightSums();

	int clauses=0;
	for (unsigned i=0; i<pi.clauses.size(); ++i) {
		if (!pi.isClauseRemoved(i)) ++clauses;
	}
	int vars=0;
	for (int i=0; i< pi.vars;++i) {
		if (!pi.isVarRemoved(i)) ++vars;
	}
	for (int i=0; i<pi.objectives; ++i) {
		stats["final_weight_sum"+std::to_string(i)]=rLog.weightRange[i];
	}
	stats["final_clauses"]=clauses;
	stats["final_variables"]=vars;
}

std::string Preprocessor::version(int l) {
	std::stringstream vs;
	if (l&1) vs << "MaxPRE ";
	vs << "2.1.0";
	if (l&2) vs << " (" << GIT_IDENTIFIER << ", " << __DATE__ << " " << __TIME__ << ")";
	return vs.str();
}

}
