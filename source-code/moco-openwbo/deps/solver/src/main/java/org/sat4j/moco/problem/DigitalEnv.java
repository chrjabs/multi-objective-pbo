/*******************************************************************************
 * SAT4J: a SATisfiability library for Java Copyright (C) 2004, 2012 Artois University and CNRS
 *
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the Eclipse Public License v1.0
 * which accompanies this distribution, and is available at
 *  http://www.eclipse.org/legal/epl-v10.html
 *
 * Alternatively, the contents of this file may be used under the terms of
 * either the GNU Lesser General Public License Version 2.1 or later (the
 * "LGPL"), in which case the provisions of the LGPL are applicable instead
 * of those above. If you wish to allow use of your version of this file only
 * under the terms of the LGPL, and not to allow others to use your version of
 * this file under the terms of the EPL, indicate your decision by deleting
 * the provisions above and replace them with the notice and other provisions
 * required by the LGPL. If you do not delete the provisions above, a recipient
 * may use your version of this file under the terms of the EPL or the LGPL.
 *
 * Contributors:
 *   CRIL - initial API and implementation
 *   João O'Neill Cortes, INESC
 *******************************************************************************/
package org.sat4j.moco.problem;

import java.util.SortedMap;
import java.util.TreeMap;
import java.util.Arrays;
import java.util.Iterator;
import java.util.Map.Entry;





/**
 * Class with the implementation of Digital environment, for
 * representation of numbers in arbitrary multi radix bases
 * @author Joao O'Neill Cortes
 */

public class DigitalEnv {

    /**
     *List of ratios. The last one will be repeated if needed.
     */
    private Integer[] ratios;

    /**
     * last base index, starting from 0.
     */

    private int basesN;

    /**
     *Default bases is simply binary
     */
    public DigitalEnv(){
	this.basesN = 1;
	this.ratios = new Integer[]{2};
}

    public DigitalEnv(Integer[] ratios){
	this.basesN = 1;
	this.setRatios(ratios);
    }

    public int getMSB(DigitalNumber digital){return digital.getDigitI(basesN); }

    public DigitalNumber toDigital(int value){

	SortedMap<Integer, Integer> digits = new TreeMap<Integer, Integer>();
	int i = 0;
	int base = 1;
	while(value != 0){
	    int ratio = getRatio(i++);
	    int digit = (value % ratio);
	    digits.put(base, digit);  
	    value-=digit;
	    value/=ratio;
	    base *= ratio;
	}
	if(i > this.basesN )
	    this.basesN = i;
	if(digits.size()==0)
	    digits.put(1, 0);
	return new DigitalNumber(digits);
    }
    public int toInt(DigitalNumber number){
	DigitalNumber.IteratorJumps iterator = number.iterator();
	int result = 0;
	while(iterator.hasNext()){
	    result += iterator.next() * iterator.currentBase();
	}
	return result;    
    }
    
    public void setRatios(Integer[] ratios){this.ratios = ratios;}
    public void setBasesN(int basesN){this.basesN = basesN;}
    public int getBasesN(){return this.basesN;}
    public int getRatio(int i){
	if(this.ratios.length <= i)
	    return this.ratios[this.ratios.length -1];
	else return this.ratios[i];
    }

    /**
     *get Base element i.
     */

    public int getBase(int index){
	int result = 1;
	for(int j = 0; j < index; j ++ )
	    result*= getRatio(j);
	return result;		
    }

    /**
     *get the index of the base, starting in 1. TODO: checky if this
     *makes sense.  If not a valid base, returns -1.
     */

    public int getBaseI(int base){
	int i = 0;
	int candidate = 1;
	while(candidate < base)
	    candidate *= this.getRatio(i++);
	if(candidate == base)
	    return i+1;
	else
	    return -1;
    }

    public class DigitalNumber implements Iterable<Integer>{

	private SortedMap<Integer, Integer> digits;

	public DigitalNumber(SortedMap<Integer, Integer> digits){
	    this.digits = digits;
	}

	public DigitalNumber plusInt(int value){
	    int result = toInt(this);
	    result += value;
	    return toDigital(result);
	}
	public int MSB(){
	    Integer key = digits.lastKey();
	    if(key == null)
		return 0;
	    return digits.get(key);
	}
	public int getDigit(int base){return this.digits.get(base);}
	public int getDigitI(int i){return this.digits.getOrDefault(getBase(i),0);}

	public IteratorJumps iterator(){
	    return new IteratorJumps();	
	}
	public IteratorContiguous iterator2(){
	    return new IteratorContiguous(false);	
	}
	public IteratorContiguous iterator3(){
	    return new IteratorContiguous(true);	
	}
	public class IteratorJumps implements Iterator<Integer>{
	    int currentBase = 1;
	    Iterator<Entry<Integer, Integer>>  current = digits.entrySet().iterator();
	    public boolean hasNext(){return current.hasNext();}
	    public Integer next(){ 
		Entry<Integer, Integer> currentEntry = current.next();
		this.currentBase = currentEntry.getKey();
		return currentEntry.getValue();
	    }
	    public Integer currentBase(){return currentBase;}
	}

	public class IteratorJumpsReversed implements Iterator<Integer>{
	    int currentBase = 1;
	    Integer[] order = new Integer[digits.size()];
	    public IteratorJumpsReversed(){
		Iterator<Entry<Integer, Integer>>  iterator = digits.entrySet().iterator();
		int i = digits.size() - 1;
		while(iterator.hasNext())
		    order[i--] = iterator.next().getKey();
	    }
	    Iterator<Integer> current = Arrays.asList(order).iterator();
	    public boolean hasNext(){return current.hasNext();}
	    public Integer next(){ 
		this.currentBase = current.next();
		return digits.get(this.currentBase);
	    }
	    public Integer currentBase(){return currentBase;}
	}

	public class IteratorContiguous implements Iterator<Integer>{
	    private int iBase = 0;
	    private int currentBase;
	    private int currentRatio;
	    boolean reverse = false;

	    public IteratorContiguous(){
		this(false);
	    }

	    public IteratorContiguous(boolean reverse){
		this.reverse = reverse;
		if(reverse)
		    iBase = getBasesN() - 1;
		currentBase = getBase(iBase);
	    }
	    public boolean hasNext(){
		if(reverse)
		    if(iBase > -1) return true; else return false;
		else
		    if( iBase < basesN) return true; else return false;}
	    public Integer next(){
		Integer result = digits.get(currentBase);
		if(reverse){
		    if(iBase > 0)
		    currentBase /= getRatio(iBase - 1);
		    iBase--;
		}else{
		    currentBase *= getRatio(iBase);
		    iBase++;
		}
		if(result == null) return 0; else return result;
		
	    }
	    public int getIBase(){
		return iBase;
	    }
	    public int currentBase(){return this.currentBase;}
	}

	/**
	 *minimal number of digits necessary to represent the corresponding value
	 */

	public int size(){
	    int result = 0;
	    Integer key = digits.lastKey();
	    if(key == null)
		return 0;
	    result = getBaseI(digits.lastKey());
	    return result;
	}
	/**
	 * returns the most significant digit of DigitalNumber, that
	 * is, the last nonzero digit.
	 */
	public int getMSB(){
	    int MSBase =  this.digits.lastKey();
	    return this.getDigit(MSBase);
	}
	/**
	 * returns the most significant base used in the
	 * representation of DigitalNumber, that is, the last base
	 * whose digit is nonzero.
	 */
	public int getMSBase(){
	    return this.digits.lastKey();
	}

	public String toString(){
	    String result = "";
	    Iterator<Integer> iterator = this.iterator2();
	    while(iterator.hasNext())
		result += iterator.next().toString();
	    return result;
}
    }
}

