package org.sat4j.moco;
import org.sat4j.moco.algorithm.UnsatSatMSU3;
import org.sat4j.moco.algorithm.algorithm;

public class UnsatSatGTETest extends algorithmTest {
    public UnsatSatGTETest(){};
    public algorithm instateAlgorithm(){
	return this.algCreator.create(1, "GTE", this.moco);
    }
}
