package org.sat4j.moco;
import org.sat4j.moco.algorithm.UnsatSat;
import org.sat4j.moco.algorithm.algorithm;

public class UnsatSatSDTest extends algorithmTest {
    public UnsatSatSDTest(){};
    public algorithm instateAlgorithm(){
	return this.algCreator.create(1, "SD", this.moco);
    }
}
