package org.sat4j.moco;
import org.sat4j.moco.algorithm.UnsatSatMSU3;
import org.sat4j.moco.algorithm.algorithm;

public class UnsatSatGTEMSU3Test extends algorithmTest {
    public UnsatSatGTEMSU3Test(){};
    public algorithm instateAlgorithm(){
	return this.algCreator.create(3, "GTE", this.moco);
    }



}
