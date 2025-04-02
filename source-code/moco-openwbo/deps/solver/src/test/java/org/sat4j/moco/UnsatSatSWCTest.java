package  org.sat4j.moco;
import org.junit.Ignore;
import org.sat4j.moco.algorithm.UnsatSat;
import org.sat4j.moco.algorithm.algorithm;

@Ignore
public class UnsatSatSWCTest extends algorithmTest {
    public UnsatSatSWCTest(){};
    public algorithm instateAlgorithm(){
	return this.algCreator.create(1, "SWC", this.moco);
    }
}
