package org.sat4j.moco;
import org.sat4j.moco.algorithm.algorithm;

public class UnsatSatSDMSU3Test extends algorithmTest {
    public UnsatSatSDMSU3Test(){};
    public algorithm instateAlgorithm(){
	return this.algCreator.create(3, "SD", this.moco);
    }



}
