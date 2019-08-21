package IMRT_Base;

/* The function solves the problem using AMPL. To do that, the function
 * first create a script to be run from the AMPL terminal. It also creates a
 * .dat file with some specific parameter values (extraXXXXX.dat). 
 * Once the problem is solved, the function return solution  * vector "x". 
 */

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.Random;
import java.util.ArrayList;
import java.util.Arrays;

/**
 *
 * @author guille
 */
public class AMPL_Solver {
    public int organs;
    public int beams;
    public int[] bmlts;
    public int totalBmlts;
    public int[] angles;
    public int[] voxels;
    public int[] aPar;
    public int[] vPar;
    public double[] wPar; //weights
    public double[] EUD0Par;
    public double[] LB;
    public double[] UB;
    public Boolean[] isTarget;
    public double epsilon;
    //public double t;
    public double[] x;
    public String jobThreadID;
    public String solver;
    public double maxIntensity;
    //public double PTV_UB;
    
    /**
     * @param args the command line arguments
     */
    //public AMPL_Solver (Organs_bkp[] o,TreatmentPlan sol, double e) throws IOException {
    public AMPL_Solver (Organs[] o,TreatmentPlan sol, double e, String jobThreadID) 
            throws IOException {
        
        
        this.organs = o.length;
        this.beams = sol.beams;
        this.bmlts=new int[this.beams];
        this.totalBmlts = sol.beamlets;
        this.angles= new int[this.beams];
        for(int i=0; i< this.beams; i++){
            this.angles[i] = sol.selAngles[i].index;
            this.bmlts[i] = sol.selAngles[i].beamlets;
        }
        this.voxels= new int[this.organs];
        this.aPar =new int[this.organs];
        this.vPar =new int[this.organs];
        this.EUD0Par =new double[this.organs];
        this.UB =new double[this.organs];
        this.LB =new double[this.organs];
        this.wPar =new double[this.organs];
        this.isTarget = new Boolean[this.organs];
        for(int i=0; i< o.length; i++){
            this.voxels[i]= o[i].totalVoxels;
            this.aPar[i] =  o[i].a;
            this.vPar[i] =  o[i].v;
            this.EUD0Par[i] =  o[i].desiredDose;
            this.UB[i] = o[i].doseUB;
            this.LB[i] = o[i].doseLB;
            this.isTarget[i]=o[i].isTarget;
            this.wPar[i] = o[i].weight;
            //this.epsilon[i] = e[i];
            //if (o[i].isTarget){
            //    this.PTV_UB =  o[i].doseUB;
            //}
        }
        this.epsilon = e;
        //this.t= o[0].doseLB;
        this.x = new double[this.totalBmlts]; 
        this.jobThreadID = jobThreadID;
        //default Solver is ipopt
        this.solver="ipopt";
        //default MaxIntensity 400
        this.maxIntensity=400;        
    }
    
    public AMPL_Solver (Organs[] o,TreatmentPlan sol, String solver, double e, String jobThreadID) 
            throws IOException {
        
        
        this.organs = o.length;
        this.beams = sol.beams;
        this.bmlts=new int[this.beams];
        this.totalBmlts = sol.beamlets;
        this.angles= new int[this.beams];
        for(int i=0; i< this.beams; i++){
            this.angles[i] = sol.selAngles[i].index;
            this.bmlts[i] = sol.selAngles[i].beamlets;
        }
        this.voxels= new int[this.organs];
        this.aPar =new int[this.organs];
        this.vPar =new int[this.organs];
        this.EUD0Par =new double[this.organs];
        this.UB =new double[this.organs];
        this.LB =new double[this.organs];
        this.wPar =new double[this.organs];
        this.isTarget = new Boolean[this.organs];
        //this.epsilon = new double [this.organs];
        for(int i=0; i< o.length; i++){
            this.voxels[i]= o[i].totalVoxels;
            this.aPar[i] =  o[i].a;
            this.vPar[i] =  o[i].v;
            this.EUD0Par[i] =  o[i].desiredDose;
            this.UB[i] = o[i].doseUB;
            this.LB[i] = o[i].doseLB;
            this.isTarget[i]=o[i].isTarget;
            this.wPar[i] = o[i].weight;
           // this.epsilon[i] = e[i];
            //if (o[i].isTarget){
            //    this.PTV_UB =  o[i].doseUB;
            //}
        }
        this.epsilon = e;
        //this.t= o[0].doseLB;
        this.x = new double[this.totalBmlts]; 
        this.jobThreadID = jobThreadID;
        this.solver=solver;
        
        
    }
    public AMPL_Solver (Organs[] o,TreatmentPlan sol, String solver,double maxIntensity, double e, String jobThreadID) 
            throws IOException {
        
        
        this.organs = o.length;
        this.beams = sol.beams;
        this.bmlts=new int[this.beams];
        this.totalBmlts = sol.beamlets;
        this.angles= new int[this.beams];
        for(int i=0; i< this.beams; i++){
            this.angles[i] = sol.selAngles[i].index;
            this.bmlts[i] = sol.selAngles[i].beamlets;
        }
        this.voxels= new int[this.organs];
        this.aPar =new int[this.organs];
        this.vPar =new int[this.organs];
        this.EUD0Par =new double[this.organs];
        this.UB =new double[this.organs];
        this.LB =new double[this.organs];
        this.wPar =new double[this.organs];
        this.isTarget = new Boolean[this.organs];
        for(int i=0; i< o.length; i++){
            this.voxels[i]= o[i].totalVoxels;
            this.aPar[i] =  o[i].a;
            this.vPar[i] =  o[i].v;
            this.EUD0Par[i] =  o[i].desiredDose;
            this.UB[i] = o[i].doseUB;
            this.LB[i] = o[i].doseLB;
            this.isTarget[i]=o[i].isTarget;
            this.wPar[i] = o[i].weight;
            //this.epsilon[i] = e[i];
            //if (o[i].isTarget){
            //    this.PTV_UB =  o[i].doseUB;
            //}
        }
        this.epsilon = e;
        //this.t= o[0].doseLB;
        this.x = new double[this.totalBmlts]; 
        this.jobThreadID = jobThreadID;
        this.solver=solver;
        if(maxIntensity == -1){
            //default MaxIntensity 400
            this.maxIntensity=400;
        }else{
            this.maxIntensity=maxIntensity;
        }
        
    }
    
    public void generateParametersFile() throws IOException{
        String parameterFile = this.jobThreadID + "extra.dat";
        //Deleting parameter file extra.txt
              
        try{
            File file = new File(parameterFile);
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(1)");
    		}
            }
    	}catch(Exception e){
            e.printStackTrace();
        }
        Random r = new Random();
        //creating the new file
        
        BufferedWriter bwParametersFile=null;
        bwParametersFile =createBufferedWriter(parameterFile);
        
        writeLine("var x := ", bwParametersFile);
        for (int i=0;i<this.totalBmlts;i++){
            int j = i+1;
            writeLine(j + " " + r.nextDouble()*50 + "\t", bwParametersFile);
        }
        writeLine(";\n", bwParametersFile);
        
        writeLine("param a := ", bwParametersFile);
        for (int i=0;i<this.organs;i++){
            int j = i+1;
            writeLine(j + " " + this.aPar[i] + "\t", bwParametersFile);
        }
        //OAR-Target
        writeLine((this.organs+1) + " 10\t", bwParametersFile);
        writeLine(";\n", bwParametersFile);
        
        for (int i=0;i<this.organs;i++){
            int j=i+1;
            writeLine("param R" + j + " := " + this.voxels[i] + ";\n", bwParametersFile);
            
            if (this.isTarget[i]){
                writeLine("param t := " + EUD0Par[i] + ";\n", bwParametersFile);
                writeLine("param OAR_targetUB := " + this.UB[i] + ";\n", bwParametersFile);
            }else{
                writeLine("param UB" + j + " := " + this.UB[i] + ";\n", bwParametersFile);
                //writeLine("param LB" + j + " := " + this.LB[i] + ";\n", bwParametersFile);
            }
        }
        
        //writeLine("param R2 := " + this.voxels[1] + ";\n", bwParametersFile);
        //writeLine("param R3 := " + this.voxels[2] + ";\n", bwParametersFile);
        writeLine("param bmlt := " + this.totalBmlts + ";\n", bwParametersFile);
        writeLine("param epsilon := " + this.epsilon + ";\n", bwParametersFile);
        
        
        
        
        bwParametersFile.close();
        
    }
    public void generateParametersFile(double[] x) throws IOException{
        //Deleting parameter file extra.txt
        String parameterFile = this.jobThreadID + "extra.dat";
        try{
            File file = new File(parameterFile);
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(2)");
    		}
            }
    	}catch(Exception e){
    	}
        Random r = new Random();
        //creating the new file
        
        BufferedWriter bwParametersFile=null;
        bwParametersFile =createBufferedWriter(parameterFile);
        
        writeLine("var x := ", bwParametersFile);
        for (int i=0;i<this.totalBmlts;i++){
            int j = i+1;
            writeLine(j + " " + x[i] + "\t", bwParametersFile);
        }
        writeLine(";\n", bwParametersFile);
        
        writeLine("param a := ", bwParametersFile);
        for (int i=0;i<this.organs;i++){
            int j = i+1;
            writeLine(j + " " + this.aPar[i] + "\t", bwParametersFile);
        }
        //OAR-Target
        writeLine((this.organs+1) + " 10\t", bwParametersFile);
        writeLine(";\n", bwParametersFile);
        
        for (int i=0;i<this.organs;i++){
            int j=i+1;
            writeLine("param R" + j + " := " + this.voxels[i] + ";\n", bwParametersFile);
            
            if (this.isTarget[i]){
                writeLine("param t := " + EUD0Par[i] + ";\n", bwParametersFile);
                writeLine("param OAR_targetUB := " + this.UB[i] + ";\n", bwParametersFile);
            }else{
                writeLine("param UB" + j + " := " + this.UB[i] + ";\n", bwParametersFile);
                //writeLine("param LB" + j + " := " + this.LB[i] + ";\n", bwParametersFile);
            }
        }
        //writeLine("param R1 := " + this.voxels[0] + ";\n", bwParametersFile);
        //writeLine("param R2 := " + this.voxels[1] + ";\n", bwParametersFile);
        //writeLine("param R3 := " + this.voxels[2] + ";\n", bwParametersFile);
        writeLine("param bmlt := " + this.totalBmlts + ";\n", bwParametersFile);
        writeLine("param epsilon := " + this.epsilon + ";\n", bwParametersFile);
        /*
        for (int i=0;i<this.organs;i++){
            int j=i+1;
            writeLine("param epsilon"+j+" := " + this.epsilon[i] + ";\n", bwParametersFile);
        }*/
        bwParametersFile.close();
    }
    
    public void generateParametersFile_logFunction() throws IOException{
        String parameterFile = this.jobThreadID + "extraLogFunction.dat";
        //Deleting parameter file extra.txt
        
        System.out.println(parameterFile);  
        try{
            File file = new File(parameterFile);
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(3)");
    		}
            }
    	}catch(Exception e){
            e.printStackTrace();
    	}
        Random r = new Random();
        //creating the new file
        
        BufferedWriter bwParametersFile=null;
        bwParametersFile =createBufferedWriter(parameterFile);
        
        writeLine("param a := ", bwParametersFile);
        for (int i=0;i<this.organs;i++){
            int j = i+1;
            writeLine(j + " " + this.aPar[i] + "\t", bwParametersFile);
        }
        //OAR-Target
        writeLine((this.organs+1) + " 10\t", bwParametersFile);//katty se cambio 10 por 50
        writeLine(";\n", bwParametersFile);
        
        writeLine(";\n", bwParametersFile);
        
        writeLine("param v := ", bwParametersFile);
        for (int i=0;i<this.organs;i++){
            int j = i+1;
            writeLine(j + " " + this.vPar[i] + "\t", bwParametersFile);
        }
        //OAR-Target
        writeLine((this.organs+1) +" 8\t", bwParametersFile);
        writeLine(";\n", bwParametersFile);
        
        writeLine("param EUD0 := ", bwParametersFile);
        for (int i=0;i<this.organs;i++){
            int j = i+1;
            writeLine(j + " " + this.EUD0Par[i] + "\t", bwParametersFile);
        }
        writeLine(";\n", bwParametersFile);
        
        for (int i=0;i<this.organs;i++){
            int j=i+1;
            writeLine("param R" + j + " := " + this.voxels[i] + ";\n", bwParametersFile);
            
            if (this.isTarget[i]){
                writeLine("param t := " + EUD0Par[i] + ";\n", bwParametersFile);
                writeLine("param OAR_targetUB := " + this.UB[i] + ";\n", bwParametersFile);
            }else{
                writeLine("param UB" + j + " := " + this.UB[i] + ";\n", bwParametersFile);
                //writeLine("param LB" + j + " := " + this.LB[i] + ";\n", bwParametersFile);
            }
        }
        writeLine("param bmlt := " + this.totalBmlts + ";\n", bwParametersFile);
        writeLine("param epsilon := " + this.epsilon + ";\n", bwParametersFile);
        /*for (int i=0;i<this.organs;i++){
            int j=i+1;
            writeLine("param epsilon"+j+" := " + this.epsilon[i] + ";\n", bwParametersFile);
        }*/
        bwParametersFile.close();
        
    }
    
    public void generateParametersFile_logFunction(double[] x) throws IOException{
        //Deleting parameter file extra.txt
        
        String parameterFile = "./" + this.jobThreadID + "extraLogFunction.dat";
        System.out.println(parameterFile);  
        try{
            File file = new File(parameterFile);
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(4)");
    		}
            }
    	}catch(Exception e){
    	}
        Random r = new Random();
        //creating the new file
        
        
        BufferedWriter bwParametersFile=null;
        bwParametersFile =createBufferedWriter(parameterFile);
        
        writeLine("var x := ", bwParametersFile);
        for (int i=0;i<this.totalBmlts;i++){
            int j = i+1;
            writeLine(j + " " + x[i] + "\t", bwParametersFile);
        }
        writeLine(";\n", bwParametersFile);
        
        writeLine("param a := ", bwParametersFile);
        for (int i=0;i<this.organs;i++){
            int j = i+1;
            writeLine(j + " " + this.aPar[i] + "\t", bwParametersFile);
        }
        //OAR-Target
        writeLine((this.organs+1) + " 10\t", bwParametersFile);
        writeLine(";\n", bwParametersFile);
        
        writeLine("param v := ", bwParametersFile);
        for (int i=0;i<this.organs;i++){
            int j = i+1;
            writeLine(j + " " + this.vPar[i] + "\t", bwParametersFile);
        }
        //OAR-Target
        writeLine((this.organs+1) + " 8\t", bwParametersFile);
        writeLine(";\n", bwParametersFile);
        
        writeLine("param EUD0 := ", bwParametersFile);
        for (int i=0;i<this.organs;i++){
            int j = i+1;
            writeLine(j + " " + this.EUD0Par[i] + "\t", bwParametersFile);
        }
        writeLine(";\n", bwParametersFile);
        
        for (int i=0;i<this.organs;i++){
            int j=i+1;
            writeLine("param R" + j + " := " + this.voxels[i] + ";\n", bwParametersFile);
            
            if (this.isTarget[i]){
                writeLine("param t := " + EUD0Par[i] + ";\n", bwParametersFile);
                writeLine("param OAR_targetUB := " + this.UB[i] + ";\n", bwParametersFile);
            }else{
                writeLine("param UB" + j + " := " + this.UB[i] + ";\n", bwParametersFile);
                //writeLine("param LB" + j + " := " + this.LB[i] + ";\n", bwParametersFile);
            }
        }
        
        //writeLine("param R1 := " + this.voxels[0] + ";\n", bwParametersFile);
        //writeLine("param R2 := " + this.voxels[1] + ";\n", bwParametersFile);
        //writeLine("param R3 := " + this.voxels[2] + ";\n", bwParametersFile);
        writeLine("param bmlt := " + this.totalBmlts + ";\n", bwParametersFile);
        writeLine("param epsilon := " + this.epsilon + ";\n", bwParametersFile);
        /*for (int i=0;i<this.organs;i++){
            int j=i+1;
            writeLine("param epsilon"+j+" := " + this.epsilon[i] + ";\n", bwParametersFile);
        }*/
        //writeLine("param t := " + this.t + ";\n", bwParametersFile);
        
        bwParametersFile.close();
        
    }
    
    public void generateParametersFile_UlogFunction() throws IOException{
        String parameterFile = this.jobThreadID + "extraULogFunction.dat";
        //Deleting parameter file extra.txt
        
        System.out.println(parameterFile);  
        try{
            File file = new File(parameterFile);
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(3)");
    		}
            }
    	}catch(Exception e){
            e.printStackTrace();
    	}
        Random r = new Random();
        //creating the new file
        
        BufferedWriter bwParametersFile=null;
        bwParametersFile =createBufferedWriter(parameterFile);
        
        writeLine("param a := ", bwParametersFile);
        for (int i=0;i<this.organs;i++){
            int j = i+1;
            writeLine(j + " " + this.aPar[i] + "\t", bwParametersFile);
        }
        //OAR-Target
        writeLine((this.organs+1) + " 10\t", bwParametersFile);//katty se cambio 10 por 50
        writeLine(";\n", bwParametersFile);
        
        writeLine(";\n", bwParametersFile);
        
        writeLine("param v := ", bwParametersFile);
        for (int i=0;i<this.organs;i++){
            int j = i+1;
            writeLine(j + " " + this.vPar[i] + "\t", bwParametersFile);
        }
        //OAR-Target
        writeLine((this.organs+1) +" 8\t", bwParametersFile);
        writeLine(";\n", bwParametersFile);
        
        writeLine("param EUD0 := ", bwParametersFile);
        for (int i=0;i<this.organs;i++){
            int j = i+1;
            writeLine(j + " " + this.EUD0Par[i] + "\t", bwParametersFile);
        }
        writeLine(";\n", bwParametersFile);
        
        for (int i=0;i<this.organs;i++){
            int j=i+1;
            writeLine("param R" + j + " := " + this.voxels[i] + ";\n", bwParametersFile);
            
            if (this.isTarget[i]){
                writeLine("param t := " + EUD0Par[i] + ";\n", bwParametersFile);
                writeLine("param OAR_targetUB := " + this.UB[i] + ";\n", bwParametersFile);
            }else{
                writeLine("param UB" + j + " := " + this.UB[i] + ";\n", bwParametersFile);
                //writeLine("param LB" + j + " := " + this.LB[i] + ";\n", bwParametersFile);
            }
        }
        writeLine("param bmlt := " + this.totalBmlts + ";\n", bwParametersFile);
        writeLine("param epsilon := " + this.epsilon + ";\n", bwParametersFile);
        /*for (int i=0;i<this.organs;i++){
            int j=i+1;
            writeLine("param epsilon"+j+" := " + this.epsilon[i] + ";\n", bwParametersFile);
        }*/
        bwParametersFile.close();
        
    }
    
    public void generateParametersFile_UlogFunction(double[] x) throws IOException{
        //Deleting parameter file extra.txt
        
        String parameterFile = "./" + this.jobThreadID + "extraULogFunction.dat";
        System.out.println(parameterFile);  
        try{
            File file = new File(parameterFile);
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(4)");
    		}
            }
    	}catch(Exception e){
    	}
        Random r = new Random();
        //creating the new file
        
        
        BufferedWriter bwParametersFile=null;
        bwParametersFile =createBufferedWriter(parameterFile);
        
        writeLine("var x := ", bwParametersFile);
        for (int i=0;i<this.totalBmlts;i++){
            int j = i+1;
            writeLine(j + " " + x[i] + "\t", bwParametersFile);
        }
        writeLine(";\n", bwParametersFile);
        
        writeLine("param a := ", bwParametersFile);
        for (int i=0;i<this.organs;i++){
            int j = i+1;
            writeLine(j + " " + this.aPar[i] + "\t", bwParametersFile);
        }
        //OAR-Target
        writeLine((this.organs+1) + " 10\t", bwParametersFile);//se cambio 50 por 10
        writeLine(";\n", bwParametersFile);
        
        writeLine("param v := ", bwParametersFile);
        for (int i=0;i<this.organs;i++){
            int j = i+1;
            writeLine(j + " " + this.vPar[i] + "\t", bwParametersFile);
        }
        //OAR-Target
        writeLine((this.organs+1) + " 8\t", bwParametersFile);
        writeLine(";\n", bwParametersFile);
        
        writeLine("param EUD0 := ", bwParametersFile);
        for (int i=0;i<this.organs;i++){
            int j = i+1;
            writeLine(j + " " + this.EUD0Par[i] + "\t", bwParametersFile);
        }
        writeLine(";\n", bwParametersFile);
        
        for (int i=0;i<this.organs;i++){
            int j=i+1;
            writeLine("param R" + j + " := " + this.voxels[i] + ";\n", bwParametersFile);
            
            if (this.isTarget[i]){
                writeLine("param t := " + EUD0Par[i] + ";\n", bwParametersFile);
                writeLine("param OAR_targetUB := " + this.UB[i] + ";\n", bwParametersFile);
            }else{
                writeLine("param UB" + j + " := " + this.UB[i] + ";\n", bwParametersFile);
                //writeLine("param LB" + j + " := " + this.LB[i] + ";\n", bwParametersFile);
            }
        }
        
        //writeLine("param R1 := " + this.voxels[0] + ";\n", bwParametersFile);
        //writeLine("param R2 := " + this.voxels[1] + ";\n", bwParametersFile);
        //writeLine("param R3 := " + this.voxels[2] + ";\n", bwParametersFile);
        writeLine("param bmlt := " + this.totalBmlts + ";\n", bwParametersFile);
        writeLine("param epsilon := " + this.epsilon + ";\n", bwParametersFile);
        /*for (int i=0;i<this.organs;i++){
            int j=i+1;
            writeLine("param epsilon"+j+" := " + this.epsilon[i] + ";\n", bwParametersFile);
        }*/
        //writeLine("param t := " + this.t + ";\n", bwParametersFile);
        
        bwParametersFile.close();
        
    }
    
    public void generateParametersFile_convexLogFunction() throws IOException{
        String parameterFile = this.jobThreadID + "extraConvexLogFunction.dat";
        //Deleting parameter file extra.txt
        
        System.out.println(parameterFile);  
        try{
            File file = new File(parameterFile);
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(3)");
    		}
            }
    	}catch(Exception e){
            e.printStackTrace();
    	}
        Random r = new Random();
        //creating the new file
        
        BufferedWriter bwParametersFile=null;
        bwParametersFile =createBufferedWriter(parameterFile);
        
        writeLine("param a := ", bwParametersFile);
        for (int i=0;i<this.organs;i++){
            int j = i+1;
            writeLine(j + " " + this.aPar[i] + "\t", bwParametersFile);
        }
        //OAR-Target
        writeLine((this.organs+1) + " 10\t", bwParametersFile);
        writeLine(";\n", bwParametersFile);
        
        //writeLine(";\n", bwParametersFile);
        
        writeLine("param v := ", bwParametersFile);
        for (int i=0;i<this.organs;i++){
            int j = i+1;
            writeLine(j + " " + this.vPar[i] + "\t", bwParametersFile);
        }
        //OAR-Target
        writeLine((this.organs+1) +" 8\t", bwParametersFile);
        writeLine(";\n", bwParametersFile);
        
        writeLine("param EUD0 := ", bwParametersFile);
        for (int i=0;i<this.organs;i++){
            int j = i+1;
            writeLine(j + " " + this.EUD0Par[i] + "\t", bwParametersFile);
        }
        writeLine(";\n", bwParametersFile);
        
        for (int i=0;i<this.organs;i++){
            int j=i+1;
            writeLine("param R" + j + " := " + this.voxels[i] + ";\n", bwParametersFile);
            
            if (this.isTarget[i]){
                writeLine("param t := " + EUD0Par[i] + ";\n", bwParametersFile);
                writeLine("param OAR_targetUB := " + this.UB[i] + ";\n", bwParametersFile);
            }else{
                writeLine("param UB" + j + " := " + this.UB[i] + ";\n", bwParametersFile);
                //writeLine("param LB" + j + " := " + this.LB[i] + ";\n", bwParametersFile);
            }
        }
        writeLine("param bmlt := " + this.totalBmlts + ";\n", bwParametersFile);
        writeLine("param epsilon := " + this.epsilon + ";\n", bwParametersFile);
        /*for (int i=0;i<this.organs;i++){
            int j=i+1;
            writeLine("param epsilon"+j+" := " + this.epsilon[i] + ";\n", bwParametersFile);
        }*/
        bwParametersFile.close();
        
    }
    
    public void generateParametersFile_convexLogFunction(double[] x) throws IOException{
        //Deleting parameter file extra.txt
        
        String parameterFile = "./" + this.jobThreadID + "extraConvexLogFunction.dat";
        System.out.println(parameterFile);  
        try{
            File file = new File(parameterFile);
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(4)");
    		}
            }
    	}catch(Exception e){
    	}
        Random r = new Random();
        //creating the new file
        
        
        BufferedWriter bwParametersFile=null;
        bwParametersFile =createBufferedWriter(parameterFile);
        
        writeLine("var x := ", bwParametersFile);
        for (int i=0;i<this.totalBmlts;i++){
            int j = i+1;
            writeLine(j + " " + x[i] + "\t", bwParametersFile);
        }
        writeLine(";\n", bwParametersFile);
        
        writeLine("param a := ", bwParametersFile);
        for (int i=0;i<this.organs;i++){
            int j = i+1;
            writeLine(j + " " + this.aPar[i] + "\t", bwParametersFile);
        }
        //OAR-Target
        writeLine((this.organs+1) + " 10\t", bwParametersFile);
        writeLine(";\n", bwParametersFile);
        
        writeLine("param v := ", bwParametersFile);
        for (int i=0;i<this.organs;i++){
            int j = i+1;
            writeLine(j + " " + this.vPar[i] + "\t", bwParametersFile);
        }
        //OAR-Target
        writeLine((this.organs+1) + " 8\t", bwParametersFile);
        writeLine(";\n", bwParametersFile);
        
        writeLine("param EUD0 := ", bwParametersFile);
        for (int i=0;i<this.organs;i++){
            int j = i+1;
            writeLine(j + " " + this.EUD0Par[i] + "\t", bwParametersFile);
        }
        writeLine(";\n", bwParametersFile);
        
        for (int i=0;i<this.organs;i++){
            int j=i+1;
            writeLine("param R" + j + " := " + this.voxels[i] + ";\n", bwParametersFile);
            
            if (this.isTarget[i]){
                writeLine("param t := " + EUD0Par[i] + ";\n", bwParametersFile);
                writeLine("param OAR_targetUB := " + this.UB[i] + ";\n", bwParametersFile);
            }else{
                writeLine("param UB" + j + " := " + this.UB[i] + ";\n", bwParametersFile);
                //writeLine("param LB" + j + " := " + this.LB[i] + ";\n", bwParametersFile);
            }
        }
        
        //writeLine("param R1 := " + this.voxels[0] + ";\n", bwParametersFile);
        //writeLine("param R2 := " + this.voxels[1] + ";\n", bwParametersFile);
        //writeLine("param R3 := " + this.voxels[2] + ";\n", bwParametersFile);
        writeLine("param bmlt := " + this.totalBmlts + ";\n", bwParametersFile);
        writeLine("param epsilon := " + this.epsilon + ";\n", bwParametersFile);
        /*for (int i=0;i<this.organs;i++){
            int j=i+1;
            writeLine("param epsilon"+j+" := " + this.epsilon[i] + ";\n", bwParametersFile);
        }*/
        //writeLine("param t := " + this.t + ";\n", bwParametersFile);
        
        bwParametersFile.close();
        
    }
    
    public void generateParametersFile_unconstrainedWeightedSum_gEUD() throws IOException{
    String parameterFile = this.jobThreadID + "extraUnconstrainedWeightedSum.dat";
        //Deleting parameter file extra.txt
        
        System.out.println(parameterFile);  
        try{
            File file = new File(parameterFile);
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(3)");
    		}
            }
    	}catch(Exception e){
            e.printStackTrace();
    	}
        Random r = new Random();
        //creating the new file
        
        BufferedWriter bwParametersFile=null;
        bwParametersFile =createBufferedWriter(parameterFile);
        
        writeLine("param a := ", bwParametersFile);
        for (int i=0;i<this.organs;i++){
            int j = i+1;
            writeLine(j + " " + this.aPar[i] + "\t", bwParametersFile);
        }
        //OAR-Target
        writeLine((this.organs+1) + " 10\t", bwParametersFile);
        writeLine(";\n", bwParametersFile);
        
        //writeLine(";\n", bwParametersFile);
        
        writeLine("param w := ", bwParametersFile);
        for (int i=0;i<this.organs;i++){
            int j = i+1;
            writeLine(j + " " + this.wPar[i] + "\t", bwParametersFile);
        }
        //OAR-Target
        writeLine((this.organs+1) +" 8\t", bwParametersFile);
        writeLine(";\n", bwParametersFile);
        
        writeLine("param EUD0 := ", bwParametersFile);
        for (int i=0;i<this.organs;i++){
            int j = i+1;
            writeLine(j + " " + this.EUD0Par[i] + "\t", bwParametersFile);
        }
        writeLine((this.organs+1) +" 76\t", bwParametersFile);//Corregir parametrizado tambien en TreatmentPlan
        writeLine(";\n", bwParametersFile);
        
        for (int i=0;i<this.organs;i++){
            int j=i+1;
            writeLine("param R" + j + " := " + this.voxels[i] + ";\n", bwParametersFile);
            
            if (this.isTarget[i]){
                writeLine("param t := " + EUD0Par[i] + ";\n", bwParametersFile);
                writeLine("param OAR_targetUB := " + this.UB[i] + ";\n", bwParametersFile);
            }else{
                writeLine("param UB" + j + " := " + this.UB[i] + ";\n", bwParametersFile);
                //writeLine("param LB" + j + " := " + this.LB[i] + ";\n", bwParametersFile);
            }
        }
        writeLine("param bmlt := " + this.totalBmlts + ";\n", bwParametersFile);
        writeLine("param epsilon := " + this.epsilon + ";\n", bwParametersFile);
        /*for (int i=0;i<this.organs;i++){
            int j=i+1;
            writeLine("param epsilon"+j+" := " + this.epsilon[i] + ";\n", bwParametersFile);
        }*/
        bwParametersFile.close();    
        
    }
    
    public void generateParametersFile_unconstrainedWeightedSum_gEUD(double[] x) throws IOException{
        //Deleting parameter file extra.txt
        
        String parameterFile = "./" + this.jobThreadID + "extraUnconstrainedWeightedSum.dat";
        System.out.println(parameterFile);  
        try{
            File file = new File(parameterFile);
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(4)");
    		}
            }
    	}catch(Exception e){
    	}
        Random r = new Random();
        //creating the new file
        
        
        BufferedWriter bwParametersFile=null;
        bwParametersFile =createBufferedWriter(parameterFile);
        
        writeLine("var x := ", bwParametersFile);
        for (int i=0;i<this.totalBmlts;i++){
            int j = i+1;
            writeLine(j + " " + x[i] + "\t", bwParametersFile);
        }
        writeLine(";\n", bwParametersFile);
        
        writeLine("param a := ", bwParametersFile);
        for (int i=0;i<this.organs;i++){
            int j = i+1;
            writeLine(j + " " + this.aPar[i] + "\t", bwParametersFile);
        }
        //OAR-Target
        writeLine((this.organs+1) + " 10\t", bwParametersFile);
        writeLine(";\n", bwParametersFile);
        
        writeLine("param w := ", bwParametersFile);
        for (int i=0;i<this.organs;i++){
            int j = i+1;
            writeLine(j + " " + this.wPar[i] + "\t", bwParametersFile);
        }
        //OAR-Target
        writeLine((this.organs+1) + " 8\t", bwParametersFile);
        writeLine(";\n", bwParametersFile);
        
        writeLine("param EUD0 := ", bwParametersFile);
        for (int i=0;i<this.organs;i++){
            int j = i+1;
            writeLine(j + " " + this.EUD0Par[i] + "\t", bwParametersFile);
        }
       // writeLine((this.organs+1) +" 82\t", bwParametersFile);//Corregir parametrizado
        writeLine(";\n", bwParametersFile);
        
        for (int i=0;i<this.organs;i++){
            int j=i+1;
            writeLine("param R" + j + " := " + this.voxels[i] + ";\n", bwParametersFile);
            
            if (this.isTarget[i]){
                writeLine("param t := " + EUD0Par[i] + ";\n", bwParametersFile);
                writeLine("param OAR_targetUB := " + this.UB[i] + ";\n", bwParametersFile);
            }else{
                writeLine("param UB" + j + " := " + this.UB[i] + ";\n", bwParametersFile);
                //writeLine("param LB" + j + " := " + this.LB[i] + ";\n", bwParametersFile);
            }
        }
        
        //writeLine("param R1 := " + this.voxels[0] + ";\n", bwParametersFile);
        //writeLine("param R2 := " + this.voxels[1] + ";\n", bwParametersFile);
        //writeLine("param R3 := " + this.voxels[2] + ";\n", bwParametersFile);
        writeLine("param bmlt := " + this.totalBmlts + ";\n", bwParametersFile);
        writeLine("param epsilon := " + this.epsilon + ";\n", bwParametersFile);
        /*for (int i=0;i<this.organs;i++){
            int j=i+1;
            writeLine("param epsilon"+j+" := " + this.epsilon[i] + ";\n", bwParametersFile);
        }*/
        //writeLine("param t := " + this.t + ";\n", bwParametersFile);
        
        bwParametersFile.close();    
    }
    
    public void generateParametersFile_adaptiveWeightedSum() throws IOException{
        String parameterFile = this.jobThreadID + "extraAdaptiveWeightedSum.dat";
        //Deleting parameter file extra.txt
        
        System.out.println(parameterFile);  
        try{
            File file = new File(parameterFile);
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(3)");
    		}
            }
    	}catch(Exception e){
            e.printStackTrace();
    	}
        Random r = new Random();
        //creating the new file
        
        BufferedWriter bwParametersFile=null;
        bwParametersFile =createBufferedWriter(parameterFile);
        
        writeLine("param a := ", bwParametersFile);
        for (int i=0;i<this.organs;i++){
            int j = i+1;
            writeLine(j + " " + this.aPar[i] + "\t", bwParametersFile);
        }
        //OAR-Target
        writeLine((this.organs+1) + " 10\t", bwParametersFile);
        writeLine(";\n", bwParametersFile);
        
        //writeLine(";\n", bwParametersFile);
        
        writeLine("param w := ", bwParametersFile);
        for (int i=0;i<this.organs;i++){
            int j = i+1;
            writeLine(j + " " + this.wPar[i] + "\t", bwParametersFile);
        }
        //OAR-Target
        writeLine((this.organs+1) +" 8\t", bwParametersFile);
        writeLine(";\n", bwParametersFile);
        
        writeLine("param EUD0 := ", bwParametersFile);
        for (int i=0;i<this.organs;i++){
            int j = i+1;
            writeLine(j + " " + this.EUD0Par[i] + "\t", bwParametersFile);
        }
        writeLine(";\n", bwParametersFile);
        
        for (int i=0;i<this.organs;i++){
            int j=i+1;
            writeLine("param R" + j + " := " + this.voxels[i] + ";\n", bwParametersFile);
            
            if (this.isTarget[i]){
                writeLine("param t := " + EUD0Par[i] + ";\n", bwParametersFile);
                writeLine("param OAR_targetUB := " + this.UB[i] + ";\n", bwParametersFile);
            }else{
                writeLine("param UB" + j + " := " + this.UB[i] + ";\n", bwParametersFile);
                //writeLine("param LB" + j + " := " + this.LB[i] + ";\n", bwParametersFile);
            }
        }
        writeLine("param bmlt := " + this.totalBmlts + ";\n", bwParametersFile);
        writeLine("param epsilon := " + this.epsilon + ";\n", bwParametersFile);
        /*for (int i=0;i<this.organs;i++){
            int j=i+1;
            writeLine("param epsilon"+j+" := " + this.epsilon[i] + ";\n", bwParametersFile);
        }*/
        bwParametersFile.close();
        
    }
    
    public void generateParametersFile_adaptiveWeightedSum(double[] x) throws IOException{
        //Deleting parameter file extra.txt
        
        String parameterFile = "./" + this.jobThreadID + "extraAdaptiveWeightedSum.dat";
        System.out.println(parameterFile);  
        try{
            File file = new File(parameterFile);
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(4)");
    		}
            }
    	}catch(Exception e){
    	}
        Random r = new Random();
        //creating the new file
        
        
        BufferedWriter bwParametersFile=null;
        bwParametersFile =createBufferedWriter(parameterFile);
        
        writeLine("var x := ", bwParametersFile);
        for (int i=0;i<this.totalBmlts;i++){
            int j = i+1;
            writeLine(j + " " + x[i] + "\t", bwParametersFile);
        }
        writeLine(";\n", bwParametersFile);
        
        writeLine("param a := ", bwParametersFile);
        for (int i=0;i<this.organs;i++){
            int j = i+1;
            writeLine(j + " " + this.aPar[i] + "\t", bwParametersFile);
        }
        //OAR-Target
        writeLine((this.organs+1) + " 10\t", bwParametersFile);
        writeLine(";\n", bwParametersFile);
        
        writeLine("param w := ", bwParametersFile);
        for (int i=0;i<this.organs;i++){
            int j = i+1;
            writeLine(j + " " + this.wPar[i] + "\t", bwParametersFile);
        }
        //OAR-Target
        writeLine((this.organs+1) + " 8\t", bwParametersFile);
        writeLine(";\n", bwParametersFile);
        
        writeLine("param EUD0 := ", bwParametersFile);
        for (int i=0;i<this.organs;i++){
            int j = i+1;
            writeLine(j + " " + this.EUD0Par[i] + "\t", bwParametersFile);
        }
        writeLine(";\n", bwParametersFile);
        
        for (int i=0;i<this.organs;i++){
            int j=i+1;
            writeLine("param R" + j + " := " + this.voxels[i] + ";\n", bwParametersFile);
            
            if (this.isTarget[i]){
                writeLine("param t := " + EUD0Par[i] + ";\n", bwParametersFile);
                writeLine("param OAR_targetUB := " + this.UB[i] + ";\n", bwParametersFile);
            }else{
                writeLine("param UB" + j + " := " + this.UB[i] + ";\n", bwParametersFile);
                //writeLine("param LB" + j + " := " + this.LB[i] + ";\n", bwParametersFile);
            }
        }
        
        //writeLine("param R1 := " + this.voxels[0] + ";\n", bwParametersFile);
        //writeLine("param R2 := " + this.voxels[1] + ";\n", bwParametersFile);
        //writeLine("param R3 := " + this.voxels[2] + ";\n", bwParametersFile);
        writeLine("param bmlt := " + this.totalBmlts + ";\n", bwParametersFile);
        writeLine("param epsilon := " + this.epsilon + ";\n", bwParametersFile);
        /*for (int i=0;i<this.organs;i++){
            int j=i+1;
            writeLine("param epsilon"+j+" := " + this.epsilon[i] + ";\n", bwParametersFile);
        }*/
        //writeLine("param t := " + this.t + ";\n", bwParametersFile);
        
        bwParametersFile.close();
        
    }
    
     public void generateParametersFile_weightedSum() throws IOException{
        String parameterFile = this.jobThreadID + "extraWeightedSum.dat";
        //Deleting parameter file extra.txt
        
        System.out.println(parameterFile);  
        try{
            File file = new File(parameterFile);
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(3)");
    		}
            }
    	}catch(Exception e){
            e.printStackTrace();
    	}
        Random r = new Random();
        //creating the new file
        
        BufferedWriter bwParametersFile=null;
        bwParametersFile =createBufferedWriter(parameterFile);
        
        writeLine("param a := ", bwParametersFile);
        for (int i=0;i<this.organs;i++){
            int j = i+1;
            writeLine(j + " " + this.aPar[i] + "\t", bwParametersFile);
        }
        //OAR-Target
        writeLine((this.organs+1) + " 10\t", bwParametersFile);
        writeLine(";\n", bwParametersFile);
        
        //writeLine(";\n", bwParametersFile);
        
        writeLine("param w := ", bwParametersFile);
        for (int i=0;i<this.organs;i++){
            int j = i+1;
            writeLine(j + " " + this.wPar[i] + "\t", bwParametersFile);
        }
        //OAR-Target
        writeLine((this.organs+1) +" 8\t", bwParametersFile);
        writeLine(";\n", bwParametersFile);
        
        writeLine("param EUD0 := ", bwParametersFile);
        for (int i=0;i<this.organs;i++){
            int j = i+1;
            writeLine(j + " " + this.EUD0Par[i] + "\t", bwParametersFile);
        }
        writeLine(";\n", bwParametersFile);
        
        for (int i=0;i<this.organs;i++){
            int j=i+1;
            writeLine("param R" + j + " := " + this.voxels[i] + ";\n", bwParametersFile);
            
            if (this.isTarget[i]){
                writeLine("param t := " + EUD0Par[i] + ";\n", bwParametersFile);
                writeLine("param OAR_targetUB := " + this.UB[i] + ";\n", bwParametersFile);
            }else{
                writeLine("param UB" + j + " := " + this.UB[i] + ";\n", bwParametersFile);
                //writeLine("param LB" + j + " := " + this.LB[i] + ";\n", bwParametersFile);
            }
        }
        writeLine("param bmlt := " + this.totalBmlts + ";\n", bwParametersFile);
        writeLine("param epsilon := " + this.epsilon + ";\n", bwParametersFile);
        /*for (int i=0;i<this.organs;i++){
            int j=i+1;
            writeLine("param epsilon"+j+" := " + this.epsilon[i] + ";\n", bwParametersFile);
        }*/
        bwParametersFile.close();
        
    }
    
    public void generateParametersFile_weightedSum(double[] x) throws IOException{
        //Deleting parameter file extra.txt
        
        String parameterFile = "./" + this.jobThreadID + "extraWeightedSum.dat";
        System.out.println(parameterFile);  
        try{
            File file = new File(parameterFile);
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(4)");
    		}
            }
    	}catch(Exception e){
    	}
        Random r = new Random();
        //creating the new file
        
        
        BufferedWriter bwParametersFile=null;
        bwParametersFile =createBufferedWriter(parameterFile);
        
        writeLine("var x := ", bwParametersFile);
        for (int i=0;i<this.totalBmlts;i++){
            int j = i+1;
            writeLine(j + " " + x[i] + "\t", bwParametersFile);
        }
        writeLine(";\n", bwParametersFile);
        
        writeLine("param a := ", bwParametersFile);
        for (int i=0;i<this.organs;i++){
            int j = i+1;
            writeLine(j + " " + this.aPar[i] + "\t", bwParametersFile);
        }
        //OAR-Target
        writeLine((this.organs+1) + " 10\t", bwParametersFile);
        writeLine(";\n", bwParametersFile);
        
        writeLine("param w := ", bwParametersFile);
        for (int i=0;i<this.organs;i++){
            int j = i+1;
            writeLine(j + " " + this.wPar[i] + "\t", bwParametersFile);
        }
        //OAR-Target
        writeLine((this.organs+1) + " 8\t", bwParametersFile);
        writeLine(";\n", bwParametersFile);
        
        writeLine("param EUD0 := ", bwParametersFile);
        for (int i=0;i<this.organs;i++){
            int j = i+1;
            writeLine(j + " " + this.EUD0Par[i] + "\t", bwParametersFile);
        }
        writeLine(";\n", bwParametersFile);
        
        for (int i=0;i<this.organs;i++){
            int j=i+1;
            writeLine("param R" + j + " := " + this.voxels[i] + ";\n", bwParametersFile);
            
            if (this.isTarget[i]){
                writeLine("param t := " + EUD0Par[i] + ";\n", bwParametersFile);
                writeLine("param OAR_targetUB := " + this.UB[i] + ";\n", bwParametersFile);
            }else{
                writeLine("param UB" + j + " := " + this.UB[i] + ";\n", bwParametersFile);
                //writeLine("param LB" + j + " := " + this.LB[i] + ";\n", bwParametersFile);
            }
        }
        
        //writeLine("param R1 := " + this.voxels[0] + ";\n", bwParametersFile);
        //writeLine("param R2 := " + this.voxels[1] + ";\n", bwParametersFile);
        //writeLine("param R3 := " + this.voxels[2] + ";\n", bwParametersFile);
        writeLine("param bmlt := " + this.totalBmlts + ";\n", bwParametersFile);
        writeLine("param epsilon := " + this.epsilon + ";\n", bwParametersFile);
        /*for (int i=0;i<this.organs;i++){
            int j=i+1;
            writeLine("param epsilon"+j+" := " + this.epsilon[i] + ";\n", bwParametersFile);
        }*/
        //writeLine("param t := " + this.t + ";\n", bwParametersFile);
        
        bwParametersFile.close();
        
    }
    
    public void getSolution() throws FileNotFoundException, IOException{
        String dir = this.jobThreadID + "currentSol.txt";
        String[] auxReader=null;
        File f = new File(dir);
        if (f.exists()) {
            try(BufferedReader fileIn = new BufferedReader(new FileReader(f))){
                String line = "";
                line=fileIn.readLine(); //avoid first line;
                line=fileIn.readLine();
                auxReader = line.split(" ");
                int j=0;
                double[] auxX = new double[this.x.length * 2];
                while (!";".equals(auxReader[0])){

                    for (String auxReader1 : auxReader) {
                        if (!"".equals(auxReader1)) {
                            auxX[j] = Double.parseDouble(auxReader1);
                            j++;
                        }
                    }
                    line=fileIn.readLine();
                    //System.out.println(line);
                    auxReader = line.split(" ");
                }
                j=0;
                for (int i=0; i<auxX.length;i++){
                    if (i%2 == 0){
                        j=(int)auxX[i];
                    }else{
                        if (j>0){
                            x[j-1] = auxX[i];
                        }else{
                            System.err.print("Archivo Existe, pero hubo un error al leer " + jobThreadID + "currentSol.txt");
                            System.err.print("auxX = ");
                            for (int k=0; k<auxX.length;k++){
                                System.err.print(auxX[k]+ " - ");
                            }
                            System.err.println();
                            System.err.print("x = ");
                            for (int k=0; k<this.x.length;k++){
                                System.err.print(this.x[k]+ " - ");
                            }
                            x[j-1] = auxX[i];
                        }
                    }       
                }
                fileIn.close();
            }
        }else{
            System.out.println("ERROR: ./" + dir + "/ file wasn't generated");
            for(int i=0; i < this.x.length;i++){
                this.x[i]=0;
            }
        }
    }

    private BufferedWriter createBufferedWriter(String route) throws IOException {
        BufferedWriter bw;
        File auxFile = new File(route);
                if (!auxFile.exists()) {
                    bw = new BufferedWriter(new FileWriter(route));
                }else{
                    bw = new BufferedWriter(new FileWriter(route,true));
                }
        return bw;     
    }
    
    public void runUnconstrainedWeightedSum_gEUD_Solver (Organs[] o) throws IOException{
        //Deleting Solution File currentSol.txt
        try{
            File file = new File(this.jobThreadID + "currentSol.txt");
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(5)");
    		}
            }
    	}catch(Exception e){
    	}
        System.out.println("Creating " + this.jobThreadID + "scriptUnconstrainedWeightedSum.sh");  
        //Creating the script file
        String scriptFile = this.jobThreadID + "scriptUnconstrainedWeightedSum.sh";
        try{
            File file = new File(scriptFile);
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(6)");
    		}
            }
    	}catch(Exception e){
    	}
        
        //creating the new file
        
        BufferedWriter bwParametersFile=null;
        bwParametersFile =createBufferedWriter(scriptFile);
        
        writeLine("model " + this.jobThreadID +"unconstrainedWeightedSumModel.mod;\n", bwParametersFile);
        for (Organs o1 : o) {
            writeLine("data "+ this.jobThreadID + "DDM_" + o1.name + ".dat;\n", bwParametersFile);
        }
        System.out.println("PASO 1");
        //writeLine("data "+ this.jobID +"DDM_RECTUM.dat;\n", bwParametersFile);
        //writeLine("data "+ this.jobID +"DDM_BLADDER.dat;\n", bwParametersFile);
        writeLine("data "+ this.jobThreadID +"extraUnconstrainedWeightedSum.dat;\n", bwParametersFile);
        writeLine("solve;\n", bwParametersFile);
        writeLine("display x > "+ this.jobThreadID +"currentSol.txt;\n", bwParametersFile);
        
        
        for (Organs o1 : o) {
            writeLine("display ((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} "
                        + "(sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j])^a["+(o1.index +1)+"]))^(1/a["+(o1.index +1)+"])" + " > " + this.jobThreadID + "gEUD_" + o1.name + ".txt;\n", bwParametersFile); 
            writeLine(" display {i in 1 .. R"+(o1.index +1)+"} (sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j]) > " + this.jobThreadID + "DVH_" + o1.name + ".txt;\n", bwParametersFile); 
            if (o1.isTarget){
                writeLine("display " + "((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} "
                        + "(sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j])^a["+(o.length +1)+"]))^(1/a["+(o.length +1)+"]) > " + this.jobThreadID + "gEUD_" + o1.name + "_UB.txt;\n", bwParametersFile); 
            }
        }
        
        
        /*for (Organs o1 : o) {
            writeLine("display gEUD_" + o1.name + " > " + this.jobThreadID + "gEUD_" + o1.name + ".txt;\n", bwParametersFile);
            writeLine("display d_" + o1.name + " > " + this.jobThreadID + "dvh_" + o1.name + ".txt;\n", bwParametersFile); 
            if (o1.isTarget){
                writeLine("display gEUD_" + o1.name + "_UB > "+"gEUD_" + o1.name + "_UB.txt;\n", bwParametersFile);
            }
        }*/
        
         bwParametersFile.close();
        
        /*CREATING THE LOGISTIC MODEL FILE FOR AMPL*/
        
        System.out.println("Creating " + this.jobThreadID + "unconstrainedWeightedSumModel.mod");  
        //Creating the script file
        scriptFile = this.jobThreadID + "unconstrainedWeightedSumModel.mod";
        try{
            File file = new File(scriptFile);
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(6)");
    		}
            }
    	}catch(Exception e){
    	}
        //creating the new file
        bwParametersFile=null;
        bwParametersFile =createBufferedWriter(scriptFile);
        
        if (null != solver)switch (solver) {
            case "ipopt":
                writeLine("option solver ipopt; \n", bwParametersFile);
                //writeLine("options ipopt_options \"linear_solver=ma57 linear_system_scaling=mc19 wantsol=8 print_level=4 tol=0.0001\"; \n", bwParametersFile);
                writeLine("options ipopt_options \"wantsol=8 print_level=0 tol=0.0001\"; \n", bwParametersFile);
                break;
            case "knitro":
                System.out.println("Entro knitro. Solver: " + solver +". Function: 1");
                writeLine("option solver knitro; \n", bwParametersFile);
                writeLine("option knitro_options \"wantsol=8 outlev=1\"; \n", bwParametersFile);
                break;
            default:
                writeLine("option solver ipopt; \n", bwParametersFile);
                //writeLine("options ipopt_options \"linear_solver=ma57 linear_system_scaling=mc19 wantsol=8 print_level=4 tol=0.0001\"; \n", bwParametersFile);
                writeLine("options ipopt_options \"wantsol=8 print_level=0 tol=0.0001\"; \n", bwParametersFile);
                
        }
       
        for (Organs o1 : o) {
            writeLine("param R" + (o1.index +1) + "; #number of voxels of " + o1.name + "\n", bwParametersFile);
        }
        writeLine("param bmlt; 	#number of beamlets \n", bwParametersFile);
        
        for (Organs o1 : o) {
            writeLine("param ddm" + o1.name + "{1 .. R" + (o1.index +1) + ", 1 .. bmlt};\n", bwParametersFile);
        }
        for (Organs o1 : o) {
            writeLine("param UB" + (o1.index +1) + ";\n", bwParametersFile);
        }
        writeLine("param t;\n" + "param epsilon;\n" + "param OAR_targetUB;\n", bwParametersFile);
        writeLine("var x {1 .. bmlt} >= 0, <=20, default 1; \n", bwParametersFile);
        for (Organs o1 : o) {
            if(o1.isTarget){
                writeLine("var gEUD_" + o1.name + "2 >=0, default 0;\n",bwParametersFile);
            }
            writeLine("var gEUD_" + o1.name + " >=0, default 0;\n",bwParametersFile);
        }
        writeLine("param a{1 .. "+ (o.length + 1) + "}; \n", bwParametersFile);
        writeLine("param w{1 .. "+ (o.length + 1) + "}; \n", bwParametersFile);
        writeLine("param EUD0{1 .. "+ (o.length + 1) + "}; \n", bwParametersFile);
        
        
        
        /*for (Organs o1 : o) {
            writeLine("par d_" + o1.name + " {i in 1 .. R"+(o1.index +1)+"} = (sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j]); \n", bwParametersFile); 
        }*/
        
        writeLine("minimize Total_Cost: ", bwParametersFile);
        for (Organs o1 : o) {
            if (o1.isTarget){
                writeLine("w["+(o1.index+1)+"]*gEUD_"+o1.name+"^2 + w["+(o.length +1)+"]*gEUD_"+o1.name+"2^2", bwParametersFile);
            }else{
                //w[1] * gEUD_PTVHD^2 + w[2] * gEUD_RECTUM^2 + w[3] * gEUD_BLADDER^2 + w[4] * gEUD_PTVHD2^2;
                writeLine("+ w["+(o1.index+1)+"] * gEUD_"+o1.name+"^2", bwParametersFile);
            }
        }
        writeLine(";\n s.t. \n", bwParametersFile);
        
        for (Organs o1 : o) {
            if (o1.isTarget){
                writeLine(o1.name+"Const	: gEUD_"+o1.name+" =  if  (\n" +
"					((1/R"+(o1.index+1)+")*(sum {i in 1..R"+(o1.index+1)+
                        "} (sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j])^a["+(o1.index+1)+"]))^(1/a["+(o1.index+1)+"]) <= EUD0["+(o1.index+1)+"]\n" +
"				) \n" +
"				then (\n" +
"					EUD0["+(o1.index+1)+"] - ((1/R"+(o1.index+1)+")*(sum {i in 1..R"+(o1.index+1)+"} (sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j])^a["+(o1.index+1)+"]))^(1/a["+(o1.index+1)+"])\n" +
"				) \n" +
"				else 0;\n", bwParametersFile);
                writeLine(o1.name+"Const2	: gEUD_"+o1.name+"2 =  if  (\n" +
"					((1/R"+(o1.index+1)+")*(sum {i in 1..R"+(o1.index+1)+
                        "} (sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j])^a["+(o.length +1)+"]))^(1/a["+(o.length +1)+"]) <= EUD0["+(o.length +1)+"]\n" +
"				) \n" +
"				then (\n" +
"					EUD0["+(o1.index+1)+"] - ((1/R"+(o1.index+1)+")*(sum {i in 1..R"+(o1.index+1)+"} (sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j])^a["+(o.length +1)+"]))^(1/a["+(o.length +1)+"])\n" +
"				) \n" +
"				else 0;\n", bwParametersFile);
              //  writeLine("constraintOAR_Target: 	((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} (sum {j in 1..bmlt} "
               //         + "x[j]*ddm"+o1.name+"[i,j])^a["+(o.length +1)+"]))^(1/a["+(o.length +1)+"]) <=OAR_targetUB;\n", bwParametersFile);
            }else{
            writeLine(o1.name+"Const	: gEUD_"+o1.name+" = if  (\n" +
"					((1/R"+(o1.index+1)+")*(sum {i in 1..R"+(o1.index+1)+"} (sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j])^a["+(o1.index+1)+"]))^(1/a["+(o1.index+1)+"]) >= EUD0["+(o1.index+1)+"]\n" +
"				) \n" +
"				then (\n" +
"					((1/R"+(o1.index+1)+")*(sum {i in 1..R"+(o1.index+1)+"} (sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j])^a["+(o1.index+1)+"]))^(1/a["+(o1.index+1)+"])-EUD0["+(o1.index+1)+"]\n" +
"				) else 0;", bwParametersFile);
                
            }
        }
        bwParametersFile.close();
       
        //Running Process
        scriptFile = this.jobThreadID + "scriptUnconstrainedWeightedSum.sh";
        Process p = new ProcessBuilder("ampl", scriptFile).start();
        InputStream is = p.getInputStream();
        InputStreamReader isr = new InputStreamReader(is);
        System.out.print("Solving unconstrainedWeightedSum Function: ");
        //System.out.print("Angles: ");
        for(int i=0; i< this.beams; i++){
            System.out.print(" " +this.angles[i]+ " -- ");
        }
        System.out.print(" // ");
        for(int i=0; i< this.EUD0Par.length; i++){
            System.out.print(this.EUD0Par[i]+ " - ");
        }
        System.out.println();
        BufferedReader br = new BufferedReader(isr);
        String line;
        while ((line = br.readLine()) != null) {
            System.out.println(line);
        }
        br.close();
    }
    
    
    public void runLogFunction_Solver (Organs[] o) throws IOException{
        //Deleting Solution File currentSol.txt
        try{
            File file = new File(this.jobThreadID + "currentSol.txt");
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(5)");
    		}
            }
    	}catch(Exception e){
    	}
        System.out.println("Creating " + this.jobThreadID + "scriptLogFunction.sh");  
        //Creating the script file
        String scriptFile = this.jobThreadID + "scriptLogFunction.sh";
        try{
            File file = new File(scriptFile);
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(6)");
    		}
            }
    	}catch(Exception e){
    	}
        
        //creating the new file
        
        BufferedWriter bwParametersFile=null;
        bwParametersFile =createBufferedWriter(scriptFile);
        
        writeLine("model " + this.jobThreadID +"logisticModel.mod;\n", bwParametersFile);
        for (Organs o1 : o) {
            writeLine("data "+ this.jobThreadID + "DDM_" + o1.name + ".dat;\n", bwParametersFile);
        }
        System.out.println("PASO 1");
        //writeLine("data "+ this.jobID +"DDM_RECTUM.dat;\n", bwParametersFile);
        //writeLine("data "+ this.jobID +"DDM_BLADDER.dat;\n", bwParametersFile);
        writeLine("data "+ this.jobThreadID +"extraLogFunction.dat;\n", bwParametersFile);
        writeLine("solve;\n", bwParametersFile);
        writeLine("display x > "+ this.jobThreadID +"currentSol.txt;\n", bwParametersFile);
        
        
        for (Organs o1 : o) {
            writeLine("display ((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} "
                        + "(sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j])^a["+(o1.index +1)+"]))^(1/a["+(o1.index +1)+"])" + " > " + this.jobThreadID + "gEUD_" + o1.name + ".txt;\n", bwParametersFile); 
            writeLine(" display {i in 1 .. R"+(o1.index +1)+"} (sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j]) > " + this.jobThreadID + "DVH_" + o1.name + ".txt;\n", bwParametersFile); 
            if (o1.isTarget){
                writeLine("display " + "((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} "
                        + "(sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j])^a["+(o.length +1)+"]))^(1/a["+(o.length +1)+"]) > " + this.jobThreadID + "gEUD_" + o1.name + "_UB.txt;\n", bwParametersFile); 
            }
        }
        
        
        /*for (Organs o1 : o) {
            writeLine("display gEUD_" + o1.name + " > " + this.jobThreadID + "gEUD_" + o1.name + ".txt;\n", bwParametersFile);
            writeLine("display d_" + o1.name + " > " + this.jobThreadID + "dvh_" + o1.name + ".txt;\n", bwParametersFile); 
            if (o1.isTarget){
                writeLine("display gEUD_" + o1.name + "_UB > "+"gEUD_" + o1.name + "_UB.txt;\n", bwParametersFile);
            }
        }*/
        
         bwParametersFile.close();
        
        /*CREATING THE LOGISTIC MODEL FILE FOR AMPL*/
        
        System.out.println("Creating " + this.jobThreadID + "logisticModel.mod");  
        //Creating the script file
        scriptFile = this.jobThreadID + "logisticModel.mod";
        try{
            File file = new File(scriptFile);
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(6)");
    		}
            }
    	}catch(Exception e){
    	}
        //creating the new file
        bwParametersFile=null;
        bwParametersFile =createBufferedWriter(scriptFile);
        
        if (null != solver)switch (solver) {
            case "ipopt":
                writeLine("option solver ipopt; \n", bwParametersFile);
                //writeLine("options ipopt_options \"linear_solver=ma57 linear_system_scaling=mc19 wantsol=8 print_level=4 tol=0.0001\"; \n", bwParametersFile);
                writeLine("options ipopt_options \"wantsol=8 print_level=0 tol=0.0001\"; \n", bwParametersFile);
                break;
            case "knitro":
                System.out.println("Entro knitro. Solver: " + solver +". Function: 2");
                writeLine("option solver knitro; \n", bwParametersFile);
                writeLine("option knitro_options \"wantsol=8 outlev=1\"; \n", bwParametersFile);
                break;
            default:
                writeLine("option solver ipopt; \n", bwParametersFile);
                //writeLine("options ipopt_options \"linear_solver=ma57 linear_system_scaling=mc19 wantsol=8 print_level=4 tol=0.0001\"; \n", bwParametersFile);
                writeLine("options ipopt_options \"wantsol=8 print_level=0 tol=0.0001\"; \n", bwParametersFile);
                
        }
        /*writeLine("option solver ipopt; \n", bwParametersFile);
        //writeLine("options ipopt_options \"linear_solver=ma57 linear_system_scaling=mc19 wantsol=8 print_level=4 tol=0.0001\"; \n", bwParametersFile);
        writeLine("options ipopt_options \"wantsol=8 print_level=4 tol=0.0001\"; \n", bwParametersFile);*/
        
        for (Organs o1 : o) {
            writeLine("param R" + (o1.index +1) + "; #number of voxels of " + o1.name + "\n", bwParametersFile);
        }
        writeLine("param bmlt; 	#number of beamlets \n", bwParametersFile);
        
        for (Organs o1 : o) {
            writeLine("param ddm" + o1.name + "{1 .. R" + (o1.index +1) + ", 1 .. bmlt};\n", bwParametersFile);
        }
        for (Organs o1 : o) {
            writeLine("param UB" + (o1.index +1) + ";\n", bwParametersFile);
        }
        writeLine("param t;\n" + "param epsilon;\n" + "param OAR_targetUB;\n", bwParametersFile);
        writeLine("var x {1 .. bmlt} >= 0, <=20, default 1; \n", bwParametersFile);//cambiar 400 a 20
        
        writeLine("param a{1 .. "+ (o.length + 1) + "}; \n", bwParametersFile);
        writeLine("param v{1 .. "+ (o.length + 1) + "}; \n", bwParametersFile);
        writeLine("param EUD0{1 .. "+ (o.length + 1) + "}; \n", bwParametersFile);
        
        
        
        /*for (Organs o1 : o) {
            writeLine("par d_" + o1.name + " {i in 1 .. R"+(o1.index +1)+"} = (sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j]); \n", bwParametersFile); 
        }*/
        
        writeLine("minimize Total_Cost: ", bwParametersFile);
        for (Organs o1 : o) {
            if (!o1.isTarget){
                writeLine("- log((1+(((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} (sum {j in 1..bmlt} "
                        + "x[j]*ddm"+o1.name+"[i,j])^a["+(o1.index +1)+"]))^(1/a["+(o1.index +1)+"])/EUD0["+(o1.index +1)+"])^v["+(o1.index +1)+"])^-1)", bwParametersFile);
            }
        }
        writeLine(";\n s.t. \n", bwParametersFile);
        
        for (Organs o1 : o) {
            if (o1.isTarget){
                writeLine("equalityTarget: 	((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} (sum {j in 1..bmlt} "
                        + "x[j]*ddm"+o1.name+"[i,j])^a["+(o1.index +1)+"]))^(1/a["+(o1.index +1)+"]) >= t; \n", bwParametersFile);
                writeLine("constraintOAR_Target: 	((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} (sum {j in 1..bmlt} "
                        + "x[j]*ddm"+o1.name+"[i,j])^a["+(o.length +1)+"]))^(1/a["+(o.length +1)+"]) <=OAR_targetUB;\n", bwParametersFile);
            }else{
                writeLine("#OAR_UB"+(o1.index +1)+": 	((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} (sum {j in 1..bmlt} "
                        + "x[j]*ddm"+o1.name+"[i,j])^a["+(o1.index +1)+"]))^(1/a["+(o1.index +1)+"]) <= UB"+(o1.index +1)+"; \n", bwParametersFile);
            }
        }
        bwParametersFile.close();
        
        System.out.println("PASO 2");
        //Running Process
        scriptFile = this.jobThreadID + "scriptLogFunction.sh";
        Process p = new ProcessBuilder("ampl", scriptFile).start();
        InputStream is = p.getInputStream();
        InputStreamReader isr = new InputStreamReader(is);
        BufferedReader br = new BufferedReader(isr);
        String line;
        System.out.print("Solving Log Function (Wu et al., 2002): ");
        //System.out.print("Angles: ");
        for(int i=0; i< this.beams; i++){
            System.out.print(" " +this.angles[i]+ " -- ");
        }
        System.out.print(" // ");
        for(int i=0; i< this.EUD0Par.length; i++){
            System.out.print(this.EUD0Par[i]+ " - ");
        }
        System.out.println();
        
        while ((line = br.readLine()) != null) {
            System.out.println(line);
        }
         br.close();
    }
    
      
    
    public void runULogFunction_Solver (Organs[] o) throws IOException{
        //Deleting Solution File currentSol.txt
        try{
            File file = new File(this.jobThreadID + "currentSol.txt");
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(5)");
    		}
            }
    	}catch(Exception e){
    	}
        System.out.println("Creating " + this.jobThreadID + "scriptULogFunction.sh");  
        //Creating the script file
        String scriptFile = this.jobThreadID + "scriptULogFunction.sh";
        try{
            File file = new File(scriptFile);
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(6)");
    		}
            }
    	}catch(Exception e){
    	}
        
        //creating the new file
        
        BufferedWriter bwParametersFile=null;
        bwParametersFile =createBufferedWriter(scriptFile);
        
        writeLine("model " + this.jobThreadID +"UlogisticModel.mod;\n", bwParametersFile);
        for (Organs o1 : o) {
            writeLine("data "+ this.jobThreadID + "DDM_" + o1.name + ".dat;\n", bwParametersFile);
        }
        System.out.println("PASO 1");
        //writeLine("data "+ this.jobID +"DDM_RECTUM.dat;\n", bwParametersFile);
        //writeLine("data "+ this.jobID +"DDM_BLADDER.dat;\n", bwParametersFile);
        writeLine("data "+ this.jobThreadID +"extraULogFunction.dat;\n", bwParametersFile);
        writeLine("solve;\n", bwParametersFile);
        writeLine("display x > "+ this.jobThreadID +"currentSol.txt;\n", bwParametersFile);
        
        
        for (Organs o1 : o) {
            writeLine("display ((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} "
                        + "(sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j])^a["+(o1.index +1)+"]))^(1/a["+(o1.index +1)+"])" + " > " + this.jobThreadID + "gEUD_" + o1.name + ".txt;\n", bwParametersFile); 
            writeLine(" display {i in 1 .. R"+(o1.index +1)+"} (sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j]) > " + this.jobThreadID + "DVH_" + o1.name + ".txt;\n", bwParametersFile); 
            if (o1.isTarget){
                writeLine("display " + "((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} "
                        + "(sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j])^a["+(o.length +1)+"]))^(1/a["+(o.length +1)+"]) > " + this.jobThreadID + "gEUD_" + o1.name + "_UB.txt;\n", bwParametersFile); 
            }
        }
       
        
        /*for (Organs o1 : o) {
            writeLine("display gEUD_" + o1.name + " > " + this.jobThreadID + "gEUD_" + o1.name + ".txt;\n", bwParametersFile);
            writeLine("display d_" + o1.name + " > " + this.jobThreadID + "dvh_" + o1.name + ".txt;\n", bwParametersFile); 
            if (o1.isTarget){
                writeLine("display gEUD_" + o1.name + "_UB > "+"gEUD_" + o1.name + "_UB.txt;\n", bwParametersFile);
            }
        }*/
        
         bwParametersFile.close();
        
        /*CREATING THE LOGISTIC MODEL FILE FOR AMPL*/
        
        System.out.println("Creating " + this.jobThreadID + "UlogisticModel.mod");  
        //Creating the script file
        scriptFile = this.jobThreadID + "UlogisticModel.mod";
        try{
            File file = new File(scriptFile);
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(6)");
    		}
            }
    	}catch(Exception e){
    	}
        //creating the new file
        bwParametersFile=null;
        bwParametersFile =createBufferedWriter(scriptFile);
        
        if (null != solver)switch (solver) {
            case "ipopt":
                writeLine("option solver ipopt; \n", bwParametersFile);
                //writeLine("options ipopt_options \"linear_solver=ma57 linear_system_scaling=mc19 wantsol=8 print_level=4 tol=0.0001\"; \n", bwParametersFile);
                writeLine("options ipopt_options \"wantsol=8 print_level=0 tol=0.0001\"; \n", bwParametersFile);
                break;
            case "knitro":
                System.out.println("Entro knitro. Solver: " + solver +". Function: 3");
                writeLine("option solver knitro; \n", bwParametersFile);
                writeLine("option knitro_options \"wantsol=8 outlev=1\"; \n", bwParametersFile);
                break;
            default:
                writeLine("option solver ipopt; \n", bwParametersFile);
                //writeLine("options ipopt_options \"linear_solver=ma57 linear_system_scaling=mc19 wantsol=8 print_level=4 tol=0.0001\"; \n", bwParametersFile);
                writeLine("options ipopt_options \"wantsol=8 print_level=0 tol=0.0001\"; \n", bwParametersFile);
                
        }
        /*writeLine("option solver ipopt; \n", bwParametersFile);
        //writeLine("options ipopt_options \"linear_solver=ma57 linear_system_scaling=mc19 wantsol=8 print_level=4 tol=0.0001\"; \n", bwParametersFile);
        writeLine("options ipopt_options \"wantsol=8 print_level=4 tol=0.0001\"; \n", bwParametersFile);*/
        
        for (Organs o1 : o) {
            writeLine("param R" + (o1.index +1) + "; #number of voxels of " + o1.name + "\n", bwParametersFile);
        }
        writeLine("param bmlt; 	#number of beamlets \n", bwParametersFile);
        
        for (Organs o1 : o) {
            writeLine("param ddm" + o1.name + "{1 .. R" + (o1.index +1) + ", 1 .. bmlt};\n", bwParametersFile);
        }
        for (Organs o1 : o) {
            writeLine("param UB" + (o1.index +1) + ";\n", bwParametersFile);
        }
        writeLine("param t;\n" + "param epsilon;\n" + "param OAR_targetUB;\n", bwParametersFile);
        writeLine("var x {1 .. bmlt} >= 0, <=20, default 1; \n", bwParametersFile);
        
        writeLine("param a{1 .. "+ (o.length + 1) + "}; \n", bwParametersFile);
        writeLine("param v{1 .. "+ (o.length + 1) + "}; \n", bwParametersFile);
        writeLine("param EUD0{1 .. "+ (o.length + 1) + "}; \n", bwParametersFile);
        
        
        
        /*for (Organs o1 : o) {
            writeLine("par d_" + o1.name + " {i in 1 .. R"+(o1.index +1)+"} = (sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j]); \n", bwParametersFile); 
        }*/
        
        writeLine("minimize Total_Cost: ", bwParametersFile);
        for (Organs o1 : o) {
            if (o1.isTarget){
                writeLine("-log ((1+(EUD0["+(o1.index +1)+"]/(((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} (sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j])^a["+(o1.index +1)+"]))^(1/a["+(o1.index +1)+"])))^v["+(o1.index +1)+"])^-1)", bwParametersFile);
            }else{
                writeLine("- log((1+(((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} (sum {j in 1..bmlt} "
                        + "x[j]*ddm"+o1.name+"[i,j])^a["+(o1.index +1)+"]))^(1/a["+(o1.index +1)+"])/EUD0["+(o1.index +1)+"])^v["+(o1.index +1)+"])^-1)", bwParametersFile);
            }
        }
        writeLine(";\n", bwParametersFile);
        bwParametersFile.close();
        
        System.out.println("PASO 2");
        //Running Process
        scriptFile = this.jobThreadID + "scriptULogFunction.sh";
        Process p = new ProcessBuilder("ampl", scriptFile).start();
        InputStream is = p.getInputStream();
        InputStreamReader isr = new InputStreamReader(is);
        BufferedReader br = new BufferedReader(isr);
        String line;
        System.out.print("Solving U Log Function (Wu et al., 2002): ");
        //System.out.print("Angles: ");
        for(int i=0; i< this.beams; i++){
            System.out.print(" " +this.angles[i]+ " -- ");
        }
        System.out.print(" // ");
        for(int i=0; i< this.EUD0Par.length; i++){
            System.out.print(this.EUD0Par[i]+ " - ");
        }
        System.out.println();
        
        while ((line = br.readLine()) != null) {
            System.out.println(line);
        }
         br.close();
    }
    
    public void runConvexLogFunction_Solver (Organs[] o) throws IOException{
        //Deleting Solution File currentSol.txt
        try{
            File file = new File(this.jobThreadID + "currentSol.txt");
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(5)");
    		}
            }
    	}catch(Exception e){
    	}
        System.out.println("Creating " + this.jobThreadID + "scriptConvexLogFunction.sh");  
        //Creating the script file
        String scriptFile = this.jobThreadID + "scriptConvexLogFunction.sh";
        try{
            File file = new File(scriptFile);
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(6)");
    		}
            }
    	}catch(Exception e){
    	}
        
        //creating the new file
        
        BufferedWriter bwParametersFile=null;
        bwParametersFile =createBufferedWriter(scriptFile);
        
        writeLine("model " + this.jobThreadID +"convexLogModel.mod;\n", bwParametersFile);
        for (Organs o1 : o) {
            writeLine("data "+ this.jobThreadID + "DDM_" + o1.name + ".dat;\n", bwParametersFile);
        }
        System.out.println("PASO 1");
        //writeLine("data "+ this.jobID +"DDM_RECTUM.dat;\n", bwParametersFile);
        //writeLine("data "+ this.jobID +"DDM_BLADDER.dat;\n", bwParametersFile);
        writeLine("data "+ this.jobThreadID +"extraConvexLogFunction.dat;\n", bwParametersFile);
        writeLine("solve;\n", bwParametersFile);
        writeLine("display x > "+ this.jobThreadID +"currentSol.txt;\n", bwParametersFile);
        
        
        for (Organs o1 : o) {
            writeLine("display ((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} "
                        + "(sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j])^a["+(o1.index +1)+"]))^(1/a["+(o1.index +1)+"])" + " > " + this.jobThreadID + "gEUD_" + o1.name + ".txt;\n", bwParametersFile); 
            writeLine(" display {i in 1 .. R"+(o1.index +1)+"} (sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j]) > " + this.jobThreadID + "DVH_" + o1.name + ".txt;\n", bwParametersFile); 
            if (o1.isTarget){
                writeLine("display " + "((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} "
                        + "(sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j])^a["+(o.length +1)+"]))^(1/a["+(o.length +1)+"]) > " + this.jobThreadID + "gEUD_" + o1.name + "_UB.txt;\n", bwParametersFile); 
            }
        }
        
        
        /*for (Organs o1 : o) {
            writeLine("display gEUD_" + o1.name + " > " + this.jobThreadID + "gEUD_" + o1.name + ".txt;\n", bwParametersFile);
            writeLine("display d_" + o1.name + " > " + this.jobThreadID + "dvh_" + o1.name + ".txt;\n", bwParametersFile); 
            if (o1.isTarget){
                writeLine("display gEUD_" + o1.name + "_UB > "+"gEUD_" + o1.name + "_UB.txt;\n", bwParametersFile);
            }
        }*/
        
         bwParametersFile.close();
        
        /*CREATING THE LOGISTIC MODEL FILE FOR AMPL*/
        
        System.out.println("Creating " + this.jobThreadID + "convexLogModel.mod");  
        //Creating the script file
        scriptFile = this.jobThreadID + "convexLogModel.mod";
        try{
            File file = new File(scriptFile);
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(6)");
    		}
            }
    	}catch(Exception e){
    	}
        //creating the new file
        bwParametersFile=null;
        bwParametersFile =createBufferedWriter(scriptFile);
        
        if (null != solver)switch (solver) {
            case "ipopt":
                writeLine("option solver ipopt; \n", bwParametersFile);
                //writeLine("options ipopt_options \"linear_solver=ma57 linear_system_scaling=mc19 wantsol=8 print_level=4 tol=0.0001\"; \n", bwParametersFile);
                writeLine("options ipopt_options \"wantsol=8 print_level=0 tol=0.0001\"; \n", bwParametersFile);
                break;
            case "knitro":
                System.out.println("Entro knitro. Solver: " + solver +". Function: 4");
                writeLine("option solver knitro; \n", bwParametersFile);
                writeLine("option knitro_options \"wantsol=8 outlev=1\"; \n", bwParametersFile);
                break;
            default:
                writeLine("option solver ipopt; \n", bwParametersFile);
                //writeLine("options ipopt_options \"linear_solver=ma57 linear_system_scaling=mc19 wantsol=8 print_level=4 tol=0.0001\"; \n", bwParametersFile);
                writeLine("options ipopt_options \"wantsol=8 print_level=0 tol=0.0001\"; \n", bwParametersFile);
                
        }

        for (Organs o1 : o) {
            writeLine("param R" + (o1.index +1) + "; #number of voxels of " + o1.name + "\n", bwParametersFile);
        }
        writeLine("param bmlt; 	#number of beamlets \n", bwParametersFile);
        
        for (Organs o1 : o) {
            writeLine("param ddm" + o1.name + "{1 .. R" + (o1.index +1) + ", 1 .. bmlt};\n", bwParametersFile);
        }
        for (Organs o1 : o) {
            writeLine("param UB" + (o1.index +1) + ";\n", bwParametersFile);
        }
        writeLine("param t;\n" + "param epsilon;\n" + "param OAR_targetUB;\n", bwParametersFile);
        writeLine("var x {1 .. bmlt} >= 0, <=20, default 1; \n", bwParametersFile);
        
        writeLine("param a{1 .. "+ (o.length + 1) + "}; \n", bwParametersFile);
        writeLine("param v{1 .. "+ (o.length + 1) + "}; \n", bwParametersFile);
        writeLine("param EUD0{1 .. "+ (o.length + 1) + "}; \n", bwParametersFile);
        
        
        
        /*for (Organs o1 : o) {
            writeLine("par d_" + o1.name + " {i in 1 .. R"+(o1.index +1)+"} = (sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j]); \n", bwParametersFile); 
        }*/
        
        writeLine("minimize Total_Cost: ", bwParametersFile);
        for (Organs o1 : o) {
            if (!o1.isTarget){
                writeLine("- log(((((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} (sum {j in 1..bmlt} "
                        + "x[j]*ddm"+o1.name+"[i,j])^a["+(o1.index +1)+"]))^(1/a["+(o1.index +1)+"])/EUD0["+(o1.index +1)+"])^v["+(o1.index +1)+"])^-1)", bwParametersFile);
            }
        }
        writeLine(";\n s.t. \n", bwParametersFile);
        
        for (Organs o1 : o) {
            if (o1.isTarget){
                writeLine("equalityTarget: 	((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} (sum {j in 1..bmlt} "
                        + "x[j]*ddm"+o1.name+"[i,j])^a["+(o1.index +1)+"]))^(1/a["+(o1.index +1)+"]) >= t; \n", bwParametersFile);
                writeLine("constraintOAR_Target: 	((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} (sum {j in 1..bmlt} "
                        + "x[j]*ddm"+o1.name+"[i,j])^a["+(o.length +1)+"]))^(1/a["+(o.length +1)+"]) <=OAR_targetUB;\n", bwParametersFile);
            }else{
                writeLine("#OAR_UB"+(o1.index +1)+": 	((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} (sum {j in 1..bmlt} "
                        + "x[j]*ddm"+o1.name+"[i,j])^a["+(o1.index +1)+"]))^(1/a["+(o1.index +1)+"]) <= UB"+(o1.index +1)+"; \n", bwParametersFile);
            }
        }
        bwParametersFile.close();
        System.out.println("PASO 2");
        //Running Process
        scriptFile = this.jobThreadID + "scriptConvexLogFunction.sh";
        Process p = new ProcessBuilder("ampl", scriptFile).start();
        InputStream is = p.getInputStream();
        InputStreamReader isr = new InputStreamReader(is);
        BufferedReader br = new BufferedReader(isr);
        String line;
        System.out.print("Solving Convex Log Function: ");
        //System.out.print("Angles: ");
        for(int i=0; i< this.beams; i++){
            System.out.print(" " +this.angles[i]+ " -- ");
        }
        System.out.print(" // ");
        for(int i=0; i< this.EUD0Par.length; i++){
            System.out.print(this.EUD0Par[i]+ " - ");
        }
        System.out.println();
        
        while ((line = br.readLine()) != null) {
            System.out.println(line);
        }
         br.close();
    }
    
    public void runWeightedSum_Solver (Organs[] o) throws IOException{
        //Deleting Solution File currentSol.txt
        try{
            File file = new File(this.jobThreadID + "currentSol.txt");
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(5)");
    		}
            }
    	}catch(Exception e){
    	}
        System.out.println("Creating " + this.jobThreadID + "scriptWeightedSum.sh");  
        //Creating the script file
        String scriptFile = this.jobThreadID + "scriptWeightedSum.sh";
        try{
            File file = new File(scriptFile);
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(6)");
    		}
            }
    	}catch(Exception e){
    	}
        
        //creating the new file
        
        BufferedWriter bwParametersFile=null;
        bwParametersFile =createBufferedWriter(scriptFile);
        
        writeLine("model " + this.jobThreadID +"weightedSumModel.mod;\n", bwParametersFile);
        for (Organs o1 : o) {
            writeLine("data "+ this.jobThreadID + "DDM_" + o1.name + ".dat;\n", bwParametersFile);
        }
        System.out.println("PASO 1");
        //writeLine("data "+ this.jobID +"DDM_RECTUM.dat;\n", bwParametersFile);
        //writeLine("data "+ this.jobID +"DDM_BLADDER.dat;\n", bwParametersFile);
        writeLine("data "+ this.jobThreadID +"extraWeightedSum.dat;\n", bwParametersFile);
        writeLine("solve;\n", bwParametersFile);
        writeLine("display x > "+ this.jobThreadID +"currentSol.txt;\n", bwParametersFile);
        
        
        for (Organs o1 : o) {
            writeLine("display ((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} "
                        + "(sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j])^a["+(o1.index +1)+"]))^(1/a["+(o1.index +1)+"])" + " > " + this.jobThreadID + "gEUD_" + o1.name + ".txt;\n", bwParametersFile); 
            writeLine(" display {i in 1 .. R"+(o1.index +1)+"} (sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j]) > " + this.jobThreadID + "DVH_" + o1.name + ".txt;\n", bwParametersFile); 
            if (o1.isTarget){
                writeLine("display " + "((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} "
                        + "(sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j])^a["+(o.length +1)+"]))^(1/a["+(o.length +1)+"]) > " + this.jobThreadID + "gEUD_" + o1.name + "_UB.txt;\n", bwParametersFile); 
            }
        }
        
        
        /*for (Organs o1 : o) {
            writeLine("display gEUD_" + o1.name + " > " + this.jobThreadID + "gEUD_" + o1.name + ".txt;\n", bwParametersFile);
            writeLine("display d_" + o1.name + " > " + this.jobThreadID + "dvh_" + o1.name + ".txt;\n", bwParametersFile); 
            if (o1.isTarget){
                writeLine("display gEUD_" + o1.name + "_UB > "+"gEUD_" + o1.name + "_UB.txt;\n", bwParametersFile);
            }
        }*/
        
         bwParametersFile.close();
        
        /*CREATING THE LOGISTIC MODEL FILE FOR AMPL*/
        
        System.out.println("Creating " + this.jobThreadID + "weightedSumModel.mod");  
        //Creating the script file
        scriptFile = this.jobThreadID + "weightedSumModel.mod";
        try{
            File file = new File(scriptFile);
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(6)");
    		}
            }
    	}catch(Exception e){
    	}
        //creating the new file
        bwParametersFile=null;
        bwParametersFile =createBufferedWriter(scriptFile);
        
        if (null != solver)switch (solver) {
            case "ipopt":
                writeLine("option solver ipopt; \n", bwParametersFile);
                //writeLine("options ipopt_options \"linear_solver=ma57 linear_system_scaling=mc19 wantsol=8 print_level=4 tol=0.0001\"; \n", bwParametersFile);
                writeLine("options ipopt_options \"wantsol=8 print_level=0 tol=0.0001\"; \n", bwParametersFile);
                break;
            case "knitro":
                System.out.println("Entro knitro. Solver: " + solver +". Function: 5");
                writeLine("option solver knitro; \n", bwParametersFile);
                writeLine("option knitro_options \"wantsol=8 outlev=1\"; \n", bwParametersFile);
                break;
            default:
                writeLine("option solver ipopt; \n", bwParametersFile);
                //writeLine("options ipopt_options \"linear_solver=ma57 linear_system_scaling=mc19 wantsol=8 print_level=4 tol=0.0001\"; \n", bwParametersFile);
                writeLine("options ipopt_options \"wantsol=8 print_level=0 tol=0.0001\"; \n", bwParametersFile);
                
        }
       
        for (Organs o1 : o) {
            writeLine("param R" + (o1.index +1) + "; #number of voxels of " + o1.name + "\n", bwParametersFile);
        }
        writeLine("param bmlt; 	#number of beamlets \n", bwParametersFile);
        
        for (Organs o1 : o) {
            writeLine("param ddm" + o1.name + "{1 .. R" + (o1.index +1) + ", 1 .. bmlt};\n", bwParametersFile);
        }
        for (Organs o1 : o) {
            writeLine("param UB" + (o1.index +1) + ";\n", bwParametersFile);
        }
        writeLine("param t;\n" + "param epsilon;\n" + "param OAR_targetUB;\n", bwParametersFile);
        writeLine("var x {1 .. bmlt} >= 0, <=20, default 1; \n", bwParametersFile);
        
        writeLine("param a{1 .. "+ (o.length + 1) + "}; \n", bwParametersFile);
        writeLine("param w{1 .. "+ (o.length + 1) + "}; \n", bwParametersFile);
        writeLine("param EUD0{1 .. "+ (o.length + 1) + "}; \n", bwParametersFile);
        
        
        
        /*for (Organs o1 : o) {
            writeLine("par d_" + o1.name + " {i in 1 .. R"+(o1.index +1)+"} = (sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j]); \n", bwParametersFile); 
        }*/
        
        writeLine("minimize Total_Cost: ", bwParametersFile);
        for (Organs o1 : o) {
            if (!o1.isTarget){
                writeLine("+ w["+(o1.index+1)+"] * ((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} (sum {j in 1..bmlt} "
                        + "x[j]*ddm"+o1.name+"[i,j])^a["+(o1.index +1)+"]))^(1/a["+(o1.index +1)+"])", bwParametersFile);
            }
        }
        writeLine(";\n s.t. \n", bwParametersFile);
        
        for (Organs o1 : o) {
            if (o1.isTarget){
                writeLine("equalityTarget: 	((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} (sum {j in 1..bmlt} "
                        + "x[j]*ddm"+o1.name+"[i,j])^a["+(o1.index +1)+"]))^(1/a["+(o1.index +1)+"]) >= t; \n", bwParametersFile);
                writeLine("constraintOAR_Target: 	((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} (sum {j in 1..bmlt} "
                        + "x[j]*ddm"+o1.name+"[i,j])^a["+(o.length +1)+"]))^(1/a["+(o.length +1)+"]) <=OAR_targetUB;\n", bwParametersFile);
            }else{
                writeLine("#OAR_UB"+(o1.index +1)+": 	((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} (sum {j in 1..bmlt} "
                        + "x[j]*ddm"+o1.name+"[i,j])^a["+(o1.index +1)+"]))^(1/a["+(o1.index +1)+"]) <= UB"+(o1.index +1)+"; \n", bwParametersFile);
            }
        }
        bwParametersFile.close();
        
        
        
        
        
        
        System.out.println("PASO 2");
        //Running Process
        scriptFile = this.jobThreadID + "scriptWeightedSum.sh";
        Process p = new ProcessBuilder("ampl", scriptFile).start();
        InputStream is = p.getInputStream();
        InputStreamReader isr = new InputStreamReader(is);
        BufferedReader br = new BufferedReader(isr);
        String line;
        System.out.print("Solving WeightedSum Function: ");
        //System.out.print("Angles: ");
        for(int i=0; i< this.beams; i++){
            System.out.print(" " +this.angles[i]+ " -- ");
        }
        System.out.print(" // ");
        for(int i=0; i< this.EUD0Par.length; i++){
            System.out.print(this.EUD0Par[i]+ " - ");
        }
        System.out.println();
        
        while ((line = br.readLine()) != null) {
            System.out.println(line);
        }
        br.close();
    }
    
    public void runVariableWeights_Solver (Organs[] o) throws IOException{
        //Deleting Solution File currentSol.txt
        try{
            File file = new File(this.jobThreadID + "currentSol.txt");
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(5)");
    		}
            }
    	}catch(Exception e){
    	}
        System.out.println("Creating " + this.jobThreadID + "scriptVariableWeights.sh with weights ( " +this.wPar[1] +" - "+this.wPar[2] + " )");  
        //Creating the script file
        String scriptFile = this.jobThreadID + "scriptVariableWeights.sh";
        try{
            File file = new File(scriptFile);
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(6)");
    		}
            }
    	}catch(Exception e){
    	}
        
        //creating the new file
        
        BufferedWriter bwParametersFile=null;
        bwParametersFile =createBufferedWriter(scriptFile);
        
        writeLine("model " + this.jobThreadID +"variableWeightsModel.mod;\n", bwParametersFile);
        for (Organs o1 : o) {
            writeLine("data "+ this.jobThreadID + "DDM_" + o1.name + ".dat;\n", bwParametersFile);
        }
        System.out.println("PASO 1");
        //writeLine("data "+ this.jobID +"DDM_RECTUM.dat;\n", bwParametersFile);
        //writeLine("data "+ this.jobID +"DDM_BLADDER.dat;\n", bwParametersFile);
        writeLine("data "+ this.jobThreadID +"extraVariableWeights.dat;\n", bwParametersFile);
        writeLine("solve;\n", bwParametersFile);
        writeLine("display x > "+ this.jobThreadID +"currentSol.txt;\n", bwParametersFile);
        
        
        for (Organs o1 : o) {
            writeLine("display ((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} "
                        + "(sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j])^a["+(o1.index +1)+"]))^(1/a["+(o1.index +1)+"])" + " > " + this.jobThreadID + "gEUD_" + o1.name + ".txt;\n", bwParametersFile); 
            writeLine(" display {i in 1 .. R"+(o1.index +1)+"} (sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j]) > " + this.jobThreadID + "DVH_" + o1.name + ".txt;\n", bwParametersFile); 
            if (o1.isTarget){
                writeLine("display " + "((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} "
                        + "(sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j])^a["+(o.length +1)+"]))^(1/a["+(o.length +1)+"]) > " + this.jobThreadID + "gEUD_" + o1.name + "_UB.txt;\n", bwParametersFile); 
            }
        }
        
        
        /*for (Organs o1 : o) {
            writeLine("display gEUD_" + o1.name + " > " + this.jobThreadID + "gEUD_" + o1.name + ".txt;\n", bwParametersFile);
            writeLine("display d_" + o1.name + " > " + this.jobThreadID + "dvh_" + o1.name + ".txt;\n", bwParametersFile); 
            if (o1.isTarget){
                writeLine("display gEUD_" + o1.name + "_UB > "+"gEUD_" + o1.name + "_UB.txt;\n", bwParametersFile);
            }
        }*/
        
         bwParametersFile.close();
        
        /*CREATING THE LOGISTIC MODEL FILE FOR AMPL*/
        
        System.out.println("Creating " + this.jobThreadID + "variableWeightsModel.mod");  
        //Creating the script file
        scriptFile = this.jobThreadID + "variableWeightsModel.mod";
        try{
            File file = new File(scriptFile);
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(6)");
    		}
            }
    	}catch(Exception e){
    	}
        //creating the new file
        bwParametersFile=null;
        bwParametersFile =createBufferedWriter(scriptFile);
        
        if (null != solver)switch (solver) {
            case "ipopt":
                writeLine("option solver ipopt; \n", bwParametersFile);
                //writeLine("options ipopt_options \"linear_solver=ma57 linear_system_scaling=mc19 wantsol=8 print_level=4 tol=0.0001\"; \n", bwParametersFile);
                writeLine("options ipopt_options \"wantsol=8 print_level=0 tol=0.0001\"; \n", bwParametersFile);
                break;
            case "knitro":
                System.out.println("Entro knitro. Solver: " + solver +". Function: 6");
                writeLine("option solver knitro; \n", bwParametersFile);
                writeLine("option knitro_options \"wantsol=8 outlev=1\"; \n", bwParametersFile);
                break;
            default:
                writeLine("option solver ipopt; \n", bwParametersFile);
                //writeLine("options ipopt_options \"linear_solver=ma57 linear_system_scaling=mc19 wantsol=8 print_level=4 tol=0.0001\"; \n", bwParametersFile);
                writeLine("options ipopt_options \"wantsol=8 print_level=0 tol=0.0001\"; \n", bwParametersFile);
                
        }
        
        for (Organs o1 : o) {
            writeLine("param R" + (o1.index +1) + "; #number of voxels of " + o1.name + "\n", bwParametersFile);
        }
        writeLine("param bmlt; 	#number of beamlets \n", bwParametersFile);
        
        for (Organs o1 : o) {
            writeLine("param ddm" + o1.name + "{1 .. R" + (o1.index +1) + ", 1 .. bmlt};\n", bwParametersFile);
        }
        for (Organs o1 : o) {
            writeLine("param UB" + (o1.index +1) + ";\n", bwParametersFile);
        }
        writeLine("param t;\n" + "param epsilon;\n" + "param OAR_targetUB;\n", bwParametersFile);
        writeLine("var x {1 .. bmlt} >= 0, <=20, default 1; \n", bwParametersFile);
        
        writeLine("param a{1 .. "+ (o.length + 1) + "}; \n", bwParametersFile);
        writeLine("param w{1 .. "+ (o.length + 1) + "}; \n", bwParametersFile);
        writeLine("param EUD0{1 .. "+ (o.length + 1) + "}; \n", bwParametersFile);
        
        
        
        /*for (Organs o1 : o) {
            writeLine("par d_" + o1.name + " {i in 1 .. R"+(o1.index +1)+"} = (sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j]); \n", bwParametersFile); 
        }*/
        
        writeLine("minimize Total_Cost: ", bwParametersFile);
        for (Organs o1 : o) {
            if (!o1.isTarget){
                writeLine("+ w["+(o1.index+1)+"] * ((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} (sum {j in 1..bmlt} "
                        + "x[j]*ddm"+o1.name+"[i,j])^a["+(o1.index +1)+"]))^(1/a["+(o1.index +1)+"])", bwParametersFile);
            }
        }
        writeLine(";\n s.t. \n", bwParametersFile);
        
        for (Organs o1 : o) {
            if (o1.isTarget){
                writeLine("equalityTarget: 	((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} (sum {j in 1..bmlt} "
                        + "x[j]*ddm"+o1.name+"[i,j])^a["+(o1.index +1)+"]))^(1/a["+(o1.index +1)+"]) >= t; \n", bwParametersFile);
                writeLine("constraintOAR_Target: 	((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} (sum {j in 1..bmlt} "
                        + "x[j]*ddm"+o1.name+"[i,j])^a["+(o.length +1)+"]))^(1/a["+(o.length +1)+"]) <=OAR_targetUB;\n", bwParametersFile);
            }else{
                writeLine("#OAR_UB"+(o1.index +1)+": 	((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} (sum {j in 1..bmlt} "
                        + "x[j]*ddm"+o1.name+"[i,j])^a["+(o1.index +1)+"]))^(1/a["+(o1.index +1)+"]) <= UB"+(o1.index +1)+"; \n", bwParametersFile);
            }
        }
        bwParametersFile.close();
         
        System.out.println("PASO 2");
        //Running Process
        scriptFile = this.jobThreadID + "scriptVariableWeights.sh";
        Process p = new ProcessBuilder("ampl", scriptFile).start();
        InputStream is = p.getInputStream();
        InputStreamReader isr = new InputStreamReader(is);
        BufferedReader br = new BufferedReader(isr);
        String line;
        System.out.print("Solving VariableWeights Function: ");
        //System.out.print("Angles: ");
        for(int i=0; i< this.beams; i++){
            System.out.print(" " +this.angles[i]+ " -- ");
        }
        System.out.print(" // ");
        for(int i=0; i< this.EUD0Par.length; i++){
            System.out.print(this.EUD0Par[i]+ " - ");
        }
        System.out.println();
        
        while ((line = br.readLine()) != null) {
            System.out.println(line);
        }
        br.close();
    }
    
    public void runAdaptiveWeightedSum_Solver (Organs[] o) throws IOException{
        //Deleting Solution File currentSol.txt
        try{
            File file = new File(this.jobThreadID + "currentSol.txt");
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(5)");
    		}
            }
    	}catch(Exception e){
    	}
        System.out.println("Creating " + this.jobThreadID + "scriptAdaptiveWeightedSum.sh with weights ( " +this.wPar[1] +" - "+this.wPar[2] + " )");  
        //Creating the script file
        String scriptFile = this.jobThreadID + "scriptAdaptiveWeightedSum.sh";
        try{
            File file = new File(scriptFile);
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(6)");
    		}
            }
    	}catch(Exception e){
    	}
        
        //creating the new file
        
        BufferedWriter bwParametersFile=null;
        bwParametersFile =createBufferedWriter(scriptFile);
        
        writeLine("model " + this.jobThreadID +"adaptiveWeightedSum.mod;\n", bwParametersFile);
        for (Organs o1 : o) {
            writeLine("data "+ this.jobThreadID + "DDM_" + o1.name + ".dat;\n", bwParametersFile);
        }
        System.out.println("PASO 1");
        //writeLine("data "+ this.jobID +"DDM_RECTUM.dat;\n", bwParametersFile);
        //writeLine("data "+ this.jobID +"DDM_BLADDER.dat;\n", bwParametersFile);
        writeLine("data "+ this.jobThreadID +"extraAdaptiveWeightedSum.dat;\n", bwParametersFile);
        writeLine("solve;\n", bwParametersFile);
        writeLine("display x > "+ this.jobThreadID +"currentSol.txt;\n", bwParametersFile);
        
        
        for (Organs o1 : o) {
            writeLine("display ((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} "
                        + "(sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j])^a["+(o1.index +1)+"]))^(1/a["+(o1.index +1)+"])" + " > " + this.jobThreadID + "gEUD_" + o1.name + ".txt;\n", bwParametersFile); 
            writeLine(" display {i in 1 .. R"+(o1.index +1)+"} (sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j]) > " + this.jobThreadID + "DVH_" + o1.name + ".txt;\n", bwParametersFile); 
            if (o1.isTarget){
                writeLine("display " + "((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} "
                        + "(sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j])^a["+(o.length +1)+"]))^(1/a["+(o.length +1)+"]) > " + this.jobThreadID + "gEUD_" + o1.name + "_UB.txt;\n", bwParametersFile); 
            }
        }
        
         bwParametersFile.close();
        
        /*CREATING THE WEIGHTED SUM MODEL FILE FOR AMPL*/
        
        System.out.println("Creating " + this.jobThreadID + "adaptiveWeightedSum.mod");  
        //Creating the script file
        scriptFile = this.jobThreadID + "adaptiveWeightedSum.mod";
        try{
            File file = new File(scriptFile);
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(6)");
    		}
            }
    	}catch(Exception e){
    	}
        //creating the new file
        bwParametersFile=null;
        bwParametersFile =createBufferedWriter(scriptFile);
		
        if (null != solver)switch (solver) {
            case "ipopt":
                writeLine("option solver ipopt; \n", bwParametersFile);
                //writeLine("options ipopt_options \"linear_solver=ma57 linear_system_scaling=mc19 wantsol=8 print_level=4 tol=0.0001\"; \n", bwParametersFile);
                writeLine("options ipopt_options \"wantsol=8 print_level=0 tol=0.0001\"; \n", bwParametersFile);
                break;
            case "knitro":
                System.out.println("Entro knitro. Solver: " + solver +". Function: 7");
                writeLine("option solver knitro; \n", bwParametersFile);
                writeLine("option knitro_options \"wantsol=8 outlev=0\"; \n", bwParametersFile);
                break;
            default:
                writeLine("option solver ipopt; \n", bwParametersFile);
                //writeLine("options ipopt_options \"linear_solver=ma57 linear_system_scaling=mc19 wantsol=8 print_level=4 tol=0.0001\"; \n", bwParametersFile);
                writeLine("options ipopt_options \"wantsol=8 print_level=0 tol=0.0001\"; \n", bwParametersFile);
                
        }
        
        for (Organs o1 : o) {
            writeLine("param R" + (o1.index +1) + "; #number of voxels of " + o1.name + "\n", bwParametersFile);
        }
        writeLine("param bmlt; 	#number of beamlets \n", bwParametersFile);
        
        for (Organs o1 : o) {
            writeLine("param ddm" + o1.name + "{1 .. R" + (o1.index +1) + ", 1 .. bmlt};\n", bwParametersFile);
        }
        for (Organs o1 : o) {
            writeLine("param UB" + (o1.index +1) + ";\n", bwParametersFile);
        }
        writeLine("param t;\n" + "param epsilon;\n" + "param OAR_targetUB;\n", bwParametersFile);
        writeLine("var x {1 .. bmlt} >= 0, <=20, default 1; \n", bwParametersFile);
        
        writeLine("param a{1 .. "+ (o.length + 1) + "}; \n", bwParametersFile);
        writeLine("param w{1 .. "+ (o.length + 1) + "}; \n", bwParametersFile);
        writeLine("param EUD0{1 .. "+ (o.length + 1) + "}; \n", bwParametersFile);
        
        
        
        /*for (Organs o1 : o) {
            writeLine("par d_" + o1.name + " {i in 1 .. R"+(o1.index +1)+"} = (sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j]); \n", bwParametersFile); 
        }*/
        
        writeLine("minimize Total_Cost: ", bwParametersFile);
        for (Organs o1 : o) {
            if (!o1.isTarget){
                writeLine("+ w["+(o1.index+1)+"] * ((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} (sum {j in 1..bmlt} "
                        + "x[j]*ddm"+o1.name+"[i,j])^a["+(o1.index +1)+"]))^(1/a["+(o1.index +1)+"])", bwParametersFile);
            }
        }
        writeLine(";\n s.t. \n", bwParametersFile);
        
        for (Organs o1 : o) {
            if (o1.isTarget){
                writeLine("equalityTarget: 	((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} (sum {j in 1..bmlt} "
                        + "x[j]*ddm"+o1.name+"[i,j])^a["+(o1.index +1)+"]))^(1/a["+(o1.index +1)+"]) >= t; \n", bwParametersFile);
                writeLine("constraintOAR_Target: 	((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} (sum {j in 1..bmlt} "
                        + "x[j]*ddm"+o1.name+"[i,j])^a["+(o.length +1)+"]))^(1/a["+(o.length +1)+"]) <=OAR_targetUB;\n", bwParametersFile);
            }else{
                writeLine("#OAR_UB"+(o1.index +1)+": 	((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} (sum {j in 1..bmlt} "
                        + "x[j]*ddm"+o1.name+"[i,j])^a["+(o1.index +1)+"]))^(1/a["+(o1.index +1)+"]) <= UB"+(o1.index +1)+"; \n", bwParametersFile);
            }
        }
        bwParametersFile.close();
        
        
        
        
        
        
        System.out.println("PASO 2");
        //Running Process
        
        scriptFile = this.jobThreadID + "scriptAdaptiveWeightedSum.sh";
        Process p = new ProcessBuilder("ampl", scriptFile).start();
        InputStream is = p.getInputStream();
        InputStreamReader isr = new InputStreamReader(is);
        BufferedReader br = new BufferedReader(isr);
        String line;
        System.out.print("Solving Adaptive WeightedSum Function: ");
        //System.out.print("Angles: ");
        for(int i=0; i< this.beams; i++){
            System.out.print(" " +this.angles[i]+ " -- ");
        }
        System.out.print(" // EUD_0: ");
        for(int i=0; i< this.EUD0Par.length; i++){
            System.out.print(this.EUD0Par[i]+ " - ");
        }
        System.out.print(" // Weights: ");
        for(int i=0; i< this.wPar.length; i++){
            System.out.print(this.wPar[i]+ " - ");
        }
        System.out.println();
        
        while ((line = br.readLine()) != null) {
            System.out.println(line);
        }
        br.close();
    }
    
    public void runLexRectum_Solver(Organs o[]) throws IOException{
        
        //Deleting Solution File currentSol.txt
        try{
            File file = new File(this.jobThreadID + "currentSol.txt");
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(7)");
    		}
            }
    	}catch(Exception e){
    	}
        
        //Create scriptLexicoRectum.sh file
        System.out.println("Creating " + this.jobThreadID + "scriptLexicoRectum.sh");  
        String scriptFile = this.jobThreadID + "scriptLexicoRectum.sh";
        try{
            File file = new File(scriptFile);
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(9.1)");
    		}
            }
    	}catch(Exception e){
    	}
        BufferedWriter bwParametersFile=null;
        bwParametersFile =createBufferedWriter(scriptFile);
        writeLine("model "+this.jobThreadID+"lexicoRectumModel.mod;\n", bwParametersFile);
        for (Organs o1 : o) {
            writeLine("data "+ this.jobThreadID + "DDM_" + o1.name + ".dat;\n", bwParametersFile);
        }
        
        writeLine("data "+ this.jobThreadID +"extra.dat;\n", bwParametersFile);
        writeLine("solve;\n", bwParametersFile);
        writeLine("display x > "+ this.jobThreadID +"currentSol.txt;\n", bwParametersFile);

        for (Organs o1 : o) {
            writeLine("display ((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} "
                        + "(sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j])^a["+(o1.index +1)+"]))^(1/a["+(o1.index +1)+"])" + " > " + this.jobThreadID + "gEUD_" + o1.name + ".txt;\n", bwParametersFile); 
            writeLine(" display {i in 1 .. R"+(o1.index +1)+"} (sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j]) > " + this.jobThreadID + "DVH_" + o1.name + ".txt;\n", bwParametersFile); 
            if (o1.isTarget){
                writeLine("display " + "((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} "
                        + "(sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j])^a["+(o.length +1)+"]))^(1/a["+(o.length +1)+"]) > " + this.jobThreadID + "gEUD_" + o1.name + "_UB.txt;\n", bwParametersFile); 
            }
        }
        writeLine("display _con > "+ this.jobThreadID +"lagrangeMultiplier.txt;\n", bwParametersFile);
        bwParametersFile.close();
        
        //Create LexicoRectum.mod for AMPL*/
        System.out.println("Creating " + this.jobThreadID + "lexicoRectumModel.mod");  
        scriptFile = this.jobThreadID + "lexicoRectumModel.mod";
        try{
            File file = new File(scriptFile);
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(9.2)");
    		}
            }
    	}catch(Exception e){
    	}
        //creating the new file
        bwParametersFile=null;
        bwParametersFile =createBufferedWriter(scriptFile);
		
        if (null != solver)switch (solver) {
            case "ipopt":
                writeLine("option solver ipopt; \n", bwParametersFile);
                //writeLine("options ipopt_options \"linear_solver=ma57 linear_system_scaling=mc19 wantsol=8 print_level=4 tol=0.0001\"; \n", bwParametersFile);
                writeLine("options ipopt_options \"wantsol=8 print_level=0 tol=0.0001\"; \n", bwParametersFile);
                break;
            case "knitro":
                System.out.println("Entro knitro. Solver: " + solver +". Function: 8");
                writeLine("option solver knitro; \n", bwParametersFile);
                writeLine("option knitro_options \"wantsol=8 outlev=o\"; \n", bwParametersFile);
                break;
            default:
                writeLine("option solver ipopt; \n", bwParametersFile);
                //writeLine("options ipopt_options \"linear_solver=ma57 linear_system_scaling=mc19 wantsol=8 print_level=4 tol=0.0001\"; \n", bwParametersFile);
                writeLine("options ipopt_options \"wantsol=8 print_level=0 tol=0.0001\"; \n", bwParametersFile);
                
        }
		
		for (Organs o1 : o) {
            writeLine("param R" + (o1.index +1) + "; #number of voxels of " + o1.name + "\n", bwParametersFile);
        }
        writeLine("param bmlt; 	#number of beamlets \n", bwParametersFile);
        
        for (Organs o1 : o) {
            writeLine("param ddm" + o1.name + "{1 .. R" + (o1.index +1) + ", 1 .. bmlt};\n", bwParametersFile);
        }
        for (Organs o1 : o) {
            writeLine("param UB" + (o1.index +1) + ";\n", bwParametersFile);
        }
        writeLine("param t;\n" + "param epsilon;\n" + "param OAR_targetUB;\n", bwParametersFile);
        writeLine("var x {1 .. bmlt} >= 0, <=20, default 1; \n", bwParametersFile);
        
        writeLine("param a{1 .. "+ (o.length + 1) + "}; \n", bwParametersFile);
        
        writeLine("minimize Total_Cost: ", bwParametersFile);
        
        for (Organs o1 : o) {
            if ("RECTUM".equals(o1.name)|| "RECTO".equals(o1.name)){
                writeLine("((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} "
                        + "(sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j])^a["+(o1.index +1)+"]))^(1/a["+(o1.index +1)+"]);\n", bwParametersFile); 
            }
        }
        writeLine(";\n s.t. \n", bwParametersFile);
        
        for (Organs o1 : o) {
            if (o1.isTarget){
                writeLine("equalityTarget: 	((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} (sum {j in 1..bmlt} "
                        + "x[j]*ddm"+o1.name+"[i,j])^a["+(o1.index +1)+"]))^(1/a["+(o1.index +1)+"]) >= t; \n", bwParametersFile);
                writeLine("constraintOAR_Target: 	((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} (sum {j in 1..bmlt} "
                        + "x[j]*ddm"+o1.name+"[i,j])^a["+(o.length +1)+"]))^(1/a["+(o.length +1)+"]) <=OAR_targetUB;\n", bwParametersFile);
            }
        }
        bwParametersFile.close();
            
        
        //Running Process
        scriptFile = this.jobThreadID + "scriptLexicoRectum.sh";
        Process p = new ProcessBuilder("ampl", scriptFile).start();
        InputStream is = p.getInputStream();
        InputStreamReader isr = new InputStreamReader(is);
        BufferedReader br = new BufferedReader(isr);
        String line;
        System.out.print("Solving LexicoRectum problem: ");
        //System.out.println("Angles: ");
        for(int i=0; i< this.beams; i++){
            System.out.print(" " +this.angles[i]+ " -- ");
        }
        System.out.print(" // ");
        for(int i=0; i< this.EUD0Par.length; i++){
            System.out.print(this.EUD0Par[i]+ " - ");
        }
        System.out.println();
        while ((line = br.readLine()) != null) {
            System.out.println(line);
        }
        br.close();
    }
    
    public void run_gEUD_Rectum_Solver (Organs[] o) throws IOException{
        //Deleting Solution File currentSol.txt
        try{
            File file = new File(this.jobThreadID + "currentSol.txt");
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(7)");
    		}
            }
    	}catch(Exception e){
    	}
        
        //Create scriptLexicoRectum.sh file
        System.out.println("Creating " + this.jobThreadID + "script_gEUD.sh");  
        String scriptFile = this.jobThreadID + "script_gEUD.sh";
        try{
            File file = new File(scriptFile);
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(11.1)");
    		}
            }
    	}catch(Exception e){
    	}
        BufferedWriter bwParametersFile=null;
        bwParametersFile =createBufferedWriter(scriptFile);
        writeLine("model "+this.jobThreadID+"gEUDConstrainedModel.mod;\n", bwParametersFile);
        for (Organs o1 : o) {
            writeLine("data "+ this.jobThreadID + "DDM_" + o1.name + ".dat;\n", bwParametersFile);
        }
        
        writeLine("data "+ this.jobThreadID +"extra.dat;\n", bwParametersFile);
        writeLine("solve;\n", bwParametersFile);
        writeLine("display x > "+ this.jobThreadID +"currentSol.txt;\n", bwParametersFile);

        for (Organs o1 : o) {
            writeLine("display ((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} "
                        + "(sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j])^a["+(o1.index +1)+"]))^(1/a["+(o1.index +1)+"])" + " > " + this.jobThreadID + "gEUD_" + o1.name + ".txt;\n", bwParametersFile); 
            writeLine(" display {i in 1 .. R"+(o1.index +1)+"} (sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j]) > " + this.jobThreadID + "DVH_" + o1.name + ".txt;\n", bwParametersFile); 
            if (o1.isTarget){
                writeLine("display " + "((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} "
                        + "(sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j])^a["+(o.length +1)+"]))^(1/a["+(o.length +1)+"]) > " + this.jobThreadID + "gEUD_" + o1.name + "_UB.txt;\n", bwParametersFile); 
            }
        }
        writeLine("display _con > "+ this.jobThreadID +"lagrangeMultiplier.txt;\n", bwParametersFile);
        bwParametersFile.close();
        
        //Create gEUDConstrainedModel.mod for AMPL*/
        System.out.println("Creating " + this.jobThreadID + "gEUDConstrainedModel.mod");  
        scriptFile = this.jobThreadID + "gEUDConstrainedModel.mod";
        try{
            File file = new File(scriptFile);
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(11.2)");
    		}
            }
    	}catch(Exception e){
    	}
        //creating the new file
        bwParametersFile=null;
        bwParametersFile =createBufferedWriter(scriptFile);
        
        if (null != solver)switch (solver) {
            case "ipopt":
                writeLine("option solver ipopt; \n", bwParametersFile);
                //writeLine("options ipopt_options \"linear_solver=ma57 linear_system_scaling=mc19 wantsol=8 print_level=4 tol=0.0001\"; \n", bwParametersFile);
                writeLine("options ipopt_options \"wantsol=8 print_level=0 tol=0.0001\"; \n", bwParametersFile);
                break;
            case "knitro":
                System.out.println("Entro knitro. Solver: " + solver +". Function: 9");
                writeLine("option solver knitro; \n", bwParametersFile);
                writeLine("option knitro_options \"wantsol=8 outlev=0\"; \n", bwParametersFile);
                break;
            default:
                writeLine("option solver ipopt; \n", bwParametersFile);
                //writeLine("options ipopt_options \"linear_solver=ma57 linear_system_scaling=mc19 wantsol=8 print_level=4 tol=0.0001\"; \n", bwParametersFile);
                writeLine("options ipopt_options \"wantsol=8 print_level=0 tol=0.0001\"; \n", bwParametersFile);
                
        }
        
        for (Organs o1 : o) {
            writeLine("param R" + (o1.index +1) + "; #number of voxels of " + o1.name + "\n", bwParametersFile);
        }
        writeLine("param bmlt; 	#number of beamlets \n", bwParametersFile);
        
        for (Organs o1 : o) {
            writeLine("param ddm" + o1.name + "{1 .. R" + (o1.index +1) + ", 1 .. bmlt};\n", bwParametersFile);
        }
        for (Organs o1 : o) {
            writeLine("param UB" + (o1.index +1) + ";\n", bwParametersFile);
        }
        writeLine("param t;\n" + "param epsilon;\n" + "param OAR_targetUB;\n", bwParametersFile);
        writeLine("var x {1 .. bmlt} >= 0, <=20, default 1; \n", bwParametersFile);
        
        writeLine("param a{1 .. "+ (o.length + 1) + "}; \n", bwParametersFile);
        
        writeLine("minimize Total_Cost: ", bwParametersFile);
        
        for (Organs o1 : o) {
            if ("RECTUM".equals(o1.name) || "RECTO".equals(o1.name)){
                writeLine("((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} "
                        + "(sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j])^a["+(o1.index +1)+"]))^(1/a["+(o1.index +1)+"]);\n", bwParametersFile); 
            }
        }
        writeLine(";\n s.t. \n", bwParametersFile);
        
        for (Organs o1 : o) {
            if (o1.isTarget){
                writeLine("equalityTarget: 	((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} (sum {j in 1..bmlt} "
                        + "x[j]*ddm"+o1.name+"[i,j])^a["+(o1.index +1)+"]))^(1/a["+(o1.index +1)+"]) >= t; \n", bwParametersFile);
                writeLine("constraintOAR_Target: 	((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} (sum {j in 1..bmlt} "
                        + "x[j]*ddm"+o1.name+"[i,j])^a["+(o.length +1)+"]))^(1/a["+(o.length +1)+"]) <=OAR_targetUB;\n", bwParametersFile);
            }else{
                if ("BLADDER".equals(o1.name) || "VEJIGA".equals(o1.name)){
                    //writeLine("constraintOAR_OAR:((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} "
                    //    + "(sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j])^a["+(o1.index +1)+"]))^(1/a["+(o1.index +1)+"]) <= epsilon;\n", bwParametersFile); 
                    writeLine("constraintOAR_OAR:((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} "
                        + "(sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j])^a["+(o1.index +1)+"]))^(1/a["+(o1.index +1)+"]) = epsilon;\n", bwParametersFile); 
                }
            }
        }
        bwParametersFile.close();
            
        
        //Running Process
        scriptFile = this.jobThreadID + "script_gEUD.sh";
        Process p = new ProcessBuilder("ampl", scriptFile).start();
        InputStream is = p.getInputStream();
        InputStreamReader isr = new InputStreamReader(is);
        BufferedReader br = new BufferedReader(isr);
        String line;
        System.out.print("Solving gEUD_Constrained problem: ");
        //System.out.println("Angles: ");
        for(int i=0; i< this.beams; i++){
            System.out.print(" " +this.angles[i]+ " -- ");
        }
        System.out.print(" // ");
        for(int i=0; i< this.EUD0Par.length; i++){
            System.out.print(this.EUD0Par[i]+ " - ");
        }
        System.out.println();
        while ((line = br.readLine()) != null) {
            System.out.println(line);
        }
        br.close();
    }
    
    public void run_gEUD_Bladder_Solver (Organs[] o) throws IOException{
        //Deleting Solution File currentSol.txt
        try{
            File file = new File(this.jobThreadID + "currentSol.txt");
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(7)");
    		}
            }
    	}catch(Exception e){
    	}
        
        //Create scriptLexicoRectum.sh file
        System.out.println("Creating " + this.jobThreadID + "script_gEUD.sh");  
        String scriptFile = this.jobThreadID + "script_gEUD.sh";
        try{
            File file = new File(scriptFile);
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(11.1)");
    		}
            }
    	}catch(Exception e){
    	}
        BufferedWriter bwParametersFile=null;
        bwParametersFile =createBufferedWriter(scriptFile);
        writeLine("model "+this.jobThreadID+"gEUDConstrainedModel.mod;\n", bwParametersFile);
        for (Organs o1 : o) {
            writeLine("data "+ this.jobThreadID + "DDM_" + o1.name + ".dat;\n", bwParametersFile);
        }
        
        writeLine("data "+ this.jobThreadID +"extra.dat;\n", bwParametersFile);
        writeLine("solve;\n", bwParametersFile);
        writeLine("display x > "+ this.jobThreadID +"currentSol.txt;\n", bwParametersFile);

        for (Organs o1 : o) {
            writeLine("display ((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} "
                        + "(sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j])^a["+(o1.index +1)+"]))^(1/a["+(o1.index +1)+"])" + " > " + this.jobThreadID + "gEUD_" + o1.name + ".txt;\n", bwParametersFile); 
            writeLine(" display {i in 1 .. R"+(o1.index +1)+"} (sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j]) > " + this.jobThreadID + "DVH_" + o1.name + ".txt;\n", bwParametersFile); 
            if (o1.isTarget){
                writeLine("display " + "((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} "
                        + "(sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j])^a["+(o.length +1)+"]))^(1/a["+(o.length +1)+"]) > " + this.jobThreadID + "gEUD_" + o1.name + "_UB.txt;\n", bwParametersFile); 
            }
        }
        writeLine("display _con > "+ this.jobThreadID +"lagrangeMultiplier.txt;\n", bwParametersFile);
        bwParametersFile.close();
        
        //Create gEUDConstrainedModel.mod for AMPL*/
        System.out.println("Creating " + this.jobThreadID + "gEUDConstrainedModel.mod");  
        scriptFile = this.jobThreadID + "gEUDConstrainedModel.mod";
        try{
            File file = new File(scriptFile);
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(11.2)");
    		}
            }
    	}catch(Exception e){
    	}
        //creating the new file
        bwParametersFile=null;
        bwParametersFile =createBufferedWriter(scriptFile);
        
        if (null != solver)switch (solver) {
            case "ipopt":
                writeLine("option solver ipopt; \n", bwParametersFile);
                //writeLine("options ipopt_options \"linear_solver=ma57 linear_system_scaling=mc19 wantsol=8 print_level=4 tol=0.0001\"; \n", bwParametersFile);
                writeLine("options ipopt_options \"wantsol=8 print_level=0 tol=0.0001\"; \n", bwParametersFile);
                break;
            case "knitro":
                System.out.println("Entro knitro. Solver: " + solver +". Function: 10");
                writeLine("option solver knitro; \n", bwParametersFile);
                writeLine("option knitro_options \"wantsol=8 outlev=1\"; \n", bwParametersFile);
                break;
            default:
                writeLine("option solver ipopt; \n", bwParametersFile);
                //writeLine("options ipopt_options \"linear_solver=ma57 linear_system_scaling=mc19 wantsol=8 print_level=4 tol=0.0001\"; \n", bwParametersFile);
                writeLine("options ipopt_options \"wantsol=8 print_level=0 tol=0.0001\"; \n", bwParametersFile);
                
        }
        
        for (Organs o1 : o) {
            writeLine("param R" + (o1.index +1) + "; #number of voxels of " + o1.name + "\n", bwParametersFile);
        }
        writeLine("param bmlt; 	#number of beamlets \n", bwParametersFile);
        
        for (Organs o1 : o) {
            writeLine("param ddm" + o1.name + "{1 .. R" + (o1.index +1) + ", 1 .. bmlt};\n", bwParametersFile);
        }
        for (Organs o1 : o) {
            writeLine("param UB" + (o1.index +1) + ";\n", bwParametersFile);
        }
        writeLine("param t;\n" + "param epsilon;\n" + "param OAR_targetUB;\n", bwParametersFile);
        writeLine("var x {1 .. bmlt} >= 0, <=20, default 1; \n", bwParametersFile);
        
        writeLine("param a{1 .. "+ (o.length + 1) + "}; \n", bwParametersFile);
        
        writeLine("minimize Total_Cost: ", bwParametersFile);
        
        for (Organs o1 : o) {
            if ("BLADDER".equals(o1.name)|| "VEJIGA".equals(o1.name)){
                writeLine("((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} "
                        + "(sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j])^a["+(o1.index +1)+"]))^(1/a["+(o1.index +1)+"]);\n", bwParametersFile); 
            }
        }
        writeLine(";\n s.t. \n", bwParametersFile);
        
        for (Organs o1 : o) {
            if (o1.isTarget){
                writeLine("equalityTarget: 	((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} (sum {j in 1..bmlt} "
                        + "x[j]*ddm"+o1.name+"[i,j])^a["+(o1.index +1)+"]))^(1/a["+(o1.index +1)+"]) >= t; \n", bwParametersFile);
                writeLine("constraintOAR_Target: 	((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} (sum {j in 1..bmlt} "
                        + "x[j]*ddm"+o1.name+"[i,j])^a["+(o.length +1)+"]))^(1/a["+(o.length +1)+"]) <=OAR_targetUB;\n", bwParametersFile);
            }else{
                if ("RECTUM".equals(o1.name) || "RECTO".equals(o1.name)){
                    //writeLine("constraintOAR_OAR:((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} "
                    //    + "(sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j])^a["+(o1.index +1)+"]))^(1/a["+(o1.index +1)+"]) <= epsilon;\n", bwParametersFile); 
                    writeLine("constraintOAR_OAR:((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} "
                        + "(sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j])^a["+(o1.index +1)+"]))^(1/a["+(o1.index +1)+"]) = epsilon;\n", bwParametersFile); 
                }
            }
        }
        bwParametersFile.close();
            
        
        //Running Process
        scriptFile = this.jobThreadID + "script_gEUD.sh";
        Process p = new ProcessBuilder("ampl", scriptFile).start();
        InputStream is = p.getInputStream();
        InputStreamReader isr = new InputStreamReader(is);
        BufferedReader br = new BufferedReader(isr);
        String line;
        System.out.print("Solving gEUD_Constrained problem: ");
        //System.out.println("Angles: ");
        for(int i=0; i< this.beams; i++){
            System.out.print(" " +this.angles[i]+ " -- ");
        }
        System.out.print(" // ");
        for(int i=0; i< this.EUD0Par.length; i++){
            System.out.print(this.EUD0Par[i]+ " - ");
        }
        System.out.println();
        while ((line = br.readLine()) != null) {
            System.out.println(line);
        }
        br.close();
    }
    
    public void runCuadraticSum_Solver (Organs[] o) throws IOException{//Maicholl
        
        //Deleting Solution File currentSol.txt
        System.out.println("Deleting Solutions Files from previous iterations for JobID "+jobThreadID+" ....Done");
        
        try{
            File file = new File("./"+this.jobThreadID + "currentSol.txt");
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(5)");
    		}
            }
            for (Organs o1 : o) {
                file = new File("./"+this.jobThreadID + "DVH_" + o1.name + ".txt");
                if (file.exists()) {
                    if(!file.delete()){
                        System.out.println("Delete operation for ./"+this.jobThreadID + "DVH_" + o1.name + ".txt failed.");
                    }
                }
                file = new File("./"+this.jobThreadID + "DVH_" + o1.name + ".txt");
                if (file.exists()) {
                    if(!file.delete()){
                        System.out.println("Delete operation for ./"+this.jobThreadID + "gEUD" + o1.name + ".txt failed.");
                    }
                }    
                if (o1.isTarget){
                    if (file.exists()) {
                        if(!file.delete()){
                            System.out.println("Delete operation for ./"+this.jobThreadID + "gEUD" + o1.name + "_UB.txt failed.");
                        }
                    }
                }
                if (file.exists()) {
                    if(!file.delete()){
                        System.out.println("Delete operation for ./"+this.jobThreadID + "voxelIndex" + o1.name + ".txt failed.");
                    }
                }
            }
    	}catch(Exception e){
            
    	}
        
        /*REPLACE AMPL FILES FOR NEXT ITERATION*/
        System.out.println("Creating " + this.jobThreadID + "scriptCuadraticSum.sh");  
        //Creating the script file
        String scriptFile = this.jobThreadID + "scriptCuadraticSum.sh";
        try{
            File file = new File(scriptFile);
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(6)");
    		}
            }
    	}catch(Exception e){
    	}
        
        //creating the new file
        
        //BufferedWriter bwParametersFile=null;
        try(BufferedWriter bwParametersFile = createBufferedWriter(scriptFile)){
        
            writeLine("model " + this.jobThreadID +"cuadraticSumModel.mod;\n", bwParametersFile);
            for (Organs o1 : o) {
                writeLine("data "+ this.jobThreadID + "DDM_" + o1.name + ".dat;\n", bwParametersFile);
            }
            writeLine("data "+ this.jobThreadID +"extraCuadraticSum.dat;\n", bwParametersFile);
            writeLine("solve;\n", bwParametersFile);
            writeLine("display x > "+ this.jobThreadID +"currentSol.txt;\n", bwParametersFile);


            for (Organs o1 : o) {
                writeLine("display ((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} "
                            + "(sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j])^a["+(o1.index +1)+"]))^(1/a["+(o1.index +1)+"])" + " > " + this.jobThreadID + "gEUD_" + o1.name + ".txt;\n", bwParametersFile); 
                writeLine(" display {i in 1 .. R"+(o1.index +1)+"} (sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j]) > " + this.jobThreadID + "DVH_" + o1.name + ".txt;\n", bwParametersFile); 
                if (o1.isTarget){
                    writeLine("display " + "((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} "
                            + "(sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j])^a["+(o.length +1)+"]))^(1/a["+(o.length +1)+"]) > " + this.jobThreadID + "gEUD_" + o1.name + "_UB.txt;\n", bwParametersFile); 
                }
            }
            bwParametersFile.close();
        }
        
        /*for (Organs o1 : o) {
            writeLine("display gEUD_" + o1.name + " > " + this.jobThreadID + "gEUD_" + o1.name + ".txt;\n", bwParametersFile);
            writeLine("display d_" + o1.name + " > " + this.jobThreadID + "dvh_" + o1.name + ".txt;\n", bwParametersFile); 
            if (o1.isTarget){
                writeLine("display gEUD_" + o1.name + "_UB > "+"gEUD_" + o1.name + "_UB.txt;\n", bwParametersFile);
            }
        }*/
        
         
        
        /*CREATING THE MODEL FILE FOR AMPL*/
        
        System.out.println("Creating " + this.jobThreadID + "cuadraticSumModel.mod");  
        //Creating the script file
        scriptFile = this.jobThreadID + "cuadraticSumModel.mod";
        try{
            File file = new File(scriptFile);
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(6)");
    		}
            }
    	}catch(Exception e){
    	}
        //creating the new file
        try(BufferedWriter bwParametersFile = createBufferedWriter(scriptFile)){
            //bwParametersFile=null;
            //bwParametersFile =createBufferedWriter(scriptFile);
           if (null != solver)switch (solver) {
                case "ipopt":
                    writeLine("option solver ipopt; \n", bwParametersFile);
                    //writeLine("options ipopt_options \"linear_solver=ma57 linear_system_scaling=mc19 wantsol=8 print_level=4 tol=0.0001\"; \n", bwParametersFile);
                    writeLine("options ipopt_options \"wantsol=8 print_level=0 tol=0.0001\"; \n", bwParametersFile);
                    break;
                case "knitro":
                    System.out.println("Entro knitro. Solver: " + solver +". Function: 11");
                    writeLine("option solver knitro; \n", bwParametersFile);
                    writeLine("option knitro_options \"wantsol=8 outlev=1\"; \n", bwParametersFile);
                    break;
                default:
                    writeLine("option solver ipopt; \n", bwParametersFile);
                    //writeLine("options ipopt_options \"linear_solver=ma57 linear_system_scaling=mc19 wantsol=8 print_level=4 tol=0.0001\"; \n", bwParametersFile);
                    writeLine("options ipopt_options \"wantsol=8 print_level=0 tol=0.0001\"; \n", bwParametersFile);
            }
        
            for (Organs o1 : o) {
                writeLine("param R" + (o1.index +1) + "; #number of voxels of " + o1.name + "\n", bwParametersFile);
            }
            writeLine("param bmlt; 	#number of beamlets \n", bwParametersFile);

            for (Organs o1 : o) {
                writeLine("param ddm" + o1.name + "{1 .. R" + (o1.index +1) + ", 1 .. bmlt};\n", bwParametersFile);
            }
            for (Organs o1 : o) {
                writeLine("#param UB" + (o1.index +1) + ";\n", bwParametersFile);
            }
            writeLine("param t;\n" + "param epsilon;\n" + "param OAR_targetUB;\n", bwParametersFile);
            writeLine("var x {1 .. bmlt} >= 0, <="+maxIntensity+", default 1; \n", bwParametersFile);

            writeLine("param a{1 .. "+ (o.length + 1) + "}; \n", bwParametersFile);
            writeLine("#param w{1 .. "+ (o.length) + "}; \n", bwParametersFile);
            writeLine("param EUD0{1 .. "+ (o.length) + "}; \n", bwParametersFile);
        
            writeLine("minimize Total_Cost: ", bwParametersFile);
            for (Organs o1 : o) {
                if (!o1.isTarget){
                   /* writeLine("+ w["+(o1.index+1)+"] * ((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} (sum {j in 1..bmlt} "
                            + "x[j]*ddm"+o1.name+"[i,j])^a["+(o1.index +1)+"]))^(1/a["+(o1.index +1)+"])", bwParametersFile);
                    */
                   writeLine("+ (1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} (max((sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j]) - EUD0["+(o1.index +1)+"],0))^2) \n", bwParametersFile); /*Maicholl*/
                }else{
                    writeLine("(1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} (max((EUD0["+(o1.index +1)+"] - (sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j])),0))^2) \n", bwParametersFile); /*Maicholl*/
                    //writeLine("+ (1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} (max((sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j]) - EUD0["+(o.length +1)+"],0)^2)) \n", bwParametersFile); /*Maicholl*/
                }
            }
       
            writeLine(";\n #s.t. \n", bwParametersFile);

            for (Organs o1 : o) {
                if (o1.isTarget){
                    writeLine("#equalityTarget: 	((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} (sum {j in 1..bmlt} "
                            + "x[j]*ddm"+o1.name+"[i,j])^a["+(o1.index +1)+"]))^(1/a["+(o1.index +1)+"]) >= t; \n", bwParametersFile);
                    writeLine("#constraintOAR_Target: 	((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} (sum {j in 1..bmlt} "
                            + "x[j]*ddm"+o1.name+"[i,j])^a["+(o.length +1)+"]))^(1/a["+(o.length +1)+"]) <=OAR_targetUB;\n", bwParametersFile);
                }else{
                    writeLine("#OAR_UB"+(o1.index +1)+": 	((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} (sum {j in 1..bmlt} "
                            + "x[j]*ddm"+o1.name+"[i,j])^a["+(o1.index +1)+"]))^(1/a["+(o1.index +1)+"]) <= UB"+(o1.index +1)+"; \n", bwParametersFile);
                }
            }
            bwParametersFile.close();
        }
        
        //Running Process
        scriptFile = this.jobThreadID + "scriptCuadraticSum.sh";
        System.out.println("Calling AMPL("+this.jobThreadID +")");
        
        Process p = new ProcessBuilder("ampl", scriptFile).start();
        InputStream is = p.getInputStream();
        InputStreamReader isr = new InputStreamReader(is);
        BufferedReader br = new BufferedReader(isr);
        String line;
        System.out.print("Solving CuadraticSum Function: ");
        //System.out.print("Angles: ");
        for(int i=0; i< this.beams; i++){
            System.out.print(" " +this.angles[i]+ " -- ");
        }
        System.out.print(" // ");
        for(int i=0; i< this.EUD0Par.length; i++){
            System.out.print(this.EUD0Par[i]+ " - ");
        }
        System.out.println();
        
        while ((line = br.readLine()) != null) {
            System.out.println(line);
        }
        is.close();
        isr.close();
        br.close();
        p.destroy();
        
    }
    
    public void runCuadraticSum_Solver (Organs[] o, int weighted) throws IOException{//Maicholl
        //Deleting Solution File currentSol.txt
        System.out.println("Deleting Solutions Files from previous iterations for JobID "+jobThreadID+" ....Done");
        int count = 0;
        try{
            File file = new File("./"+this.jobThreadID + "currentSol.txt");
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(5)");
    		}
            }
            for (Organs o1 : o) {
                file = new File("./"+this.jobThreadID + "DVH_" + o1.name + ".txt");
                if (file.exists()) {
                    if(!file.delete()){
                        System.out.println("Delete operation for ./"+this.jobThreadID + "DVH_" + o1.name + ".txt failed.");
                    }
                }
                file = new File("./"+this.jobThreadID + "DVH_" + o1.name + ".txt");
                if (file.exists()) {
                    if(!file.delete()){
                        System.out.println("Delete operation for ./"+this.jobThreadID + "gEUD" + o1.name + ".txt failed.");
                    }
                }    
                if (o1.isTarget){
                    if (file.exists()) {
                        if(!file.delete()){
                            System.out.println("Delete operation for ./"+this.jobThreadID + "gEUD" + o1.name + "_UB.txt failed.");
                        }
                    }
                }
                if (file.exists()) {
                    if(!file.delete()){
                        System.out.println("Delete operation for ./"+this.jobThreadID + "voxelIndex" + o1.name + ".txt failed.");
                    }
                }
            }
    	}catch(Exception e){
            
    	}
        
        /*REPLACE AMPL FILES FOR NEXT ITERATION*/
        System.out.println("Creating " + this.jobThreadID + "scriptCuadraticSum.sh");  
        //Creating the script file
        String scriptFile ="./"+ this.jobThreadID + "scriptCuadraticSum.sh";
        try{
            File file = new File(scriptFile);
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(6)");
    		}
            }
    	}catch(Exception e){
            e.printStackTrace();
    	}
        
        //creating the new file
        
        //BufferedWriter bwParametersFile=null;
        try(BufferedWriter bwParametersFile =createBufferedWriter(scriptFile)){
        
            writeLine("model " + this.jobThreadID +"cuadraticSumModel.mod;\n", bwParametersFile);
            for (Organs o1 : o) {
                writeLine("data "+ this.jobThreadID + "DDM_" + o1.name + ".dat;\n", bwParametersFile);
            }
            writeLine("data "+ this.jobThreadID +"extraCuadraticSum.dat;\n", bwParametersFile);
            writeLine("solve;\n", bwParametersFile);
            writeLine("display x > "+ this.jobThreadID +"currentSol.txt;\n", bwParametersFile);


            for (Organs o1 : o) {
                writeLine("display ((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} "
                            + "(sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j])^a["+(o1.index +1)+"]))^(1/a["+(o1.index +1)+"])" + " > " + this.jobThreadID + "gEUD_" + o1.name + ".txt;\n", bwParametersFile); 
                writeLine(" display {i in 1 .. R"+(o1.index +1)+"} (sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j]) > " + this.jobThreadID + "DVH_" + o1.name + ".txt;\n", bwParametersFile); 
                if (o1.isTarget){
                    writeLine("display " + "((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} "
                            + "(sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j])^a["+(o.length +1)+"]))^(1/a["+(o.length +1)+"]) > " + this.jobThreadID + "gEUD_" + o1.name + "_UB.txt;\n", bwParametersFile); 
                }
            }
            bwParametersFile.close();
        }
        /*CREATING THE LOGISTIC MODEL FILE FOR AMPL*/
        
        System.out.println("Creating " + this.jobThreadID + "cuadraticSumModel.mod");  
        //Creating the script file
        scriptFile ="./"+ this.jobThreadID + "cuadraticSumModel.mod";
        try{
            File file = new File(scriptFile);
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(6)");
    		}
            }
    	}catch(Exception e){
    	}
        //creating the new file
        //bwParametersFile=null;
        try(BufferedWriter bwParametersFile =createBufferedWriter(scriptFile)){
        
            if (null != solver)switch (solver) {
                case "ipopt":
                    writeLine("option solver ipopt; \n", bwParametersFile);
                    //writeLine("options ipopt_options \"linear_solver=ma57 linear_system_scaling=mc19 wantsol=8 print_level=4 tol=0.0001\"; \n", bwParametersFile);
                    writeLine("options ipopt_options \"wantsol=8 print_level=0 tol=0.0001\"; \n", bwParametersFile);
                    break;
                case "knitro":
                    System.out.println("Entro knitro. Solver: " + solver +". Function: 12");
                    writeLine("option solver knitro; \n", bwParametersFile);
                    writeLine("option knitro_options \"wantsol=8 outlev=1\"; \n", bwParametersFile);
                    break;
                default:
                    writeLine("option solver ipopt; \n", bwParametersFile);
                    //writeLine("options ipopt_options \"linear_solver=ma57 linear_system_scaling=mc19 wantsol=8 print_level=4 tol=0.0001\"; \n", bwParametersFile);
                    writeLine("options ipopt_options \"wantsol=8 print_level=0 tol=0.0001\"; \n", bwParametersFile);

            }

            for (Organs o1 : o) {
                writeLine("param R" + (o1.index +1) + "; #number of voxels of " + o1.name + "\n", bwParametersFile);
            }
            writeLine("param bmlt; 	#number of beamlets \n", bwParametersFile);

            for (Organs o1 : o) {
                writeLine("param ddm" + o1.name + "{1 .. R" + (o1.index +1) + ", 1 .. bmlt};\n", bwParametersFile);
            }
            for (Organs o1 : o) {
                writeLine("param UB" + (o1.index +1) + ";\n", bwParametersFile);
            }
            writeLine("param t;\n" + "param epsilon;\n" + "param OAR_targetUB;\n", bwParametersFile);
            writeLine("var x {1 .. bmlt} >= 0, <="+maxIntensity+", default 1; \n", bwParametersFile);

            writeLine("param a{1 .. "+ (o.length + 1) + "}; \n", bwParametersFile);
            writeLine("#param w{1 .. "+ (o.length) + "}; \n", bwParametersFile);
            writeLine("param EUD0{1 .. "+ (o.length) + "}; \n", bwParametersFile);



            /*for (Organs o1 : o) {
                writeLine("par d_" + o1.name + " {i in 1 .. R"+(o1.index +1)+"} = (sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j]); \n", bwParametersFile); 
            }*/

            writeLine("minimize Total_Cost: ", bwParametersFile);
            for (Organs o1 : o) {
                if (!o1.isTarget){
                    if(weighted == 0){
                        writeLine("+ (1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} (max((sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j]) - EUD0["+(o1.index +1)+"],0))^2) \n", bwParametersFile); /*Maicholl*/
                    }else{
                        writeLine("+ "+wPar[count]+"*(1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} (max((sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j]) - EUD0["+(o1.index +1)+"],0))^2) \n", bwParametersFile); /*Maicholl*/  
                    }
                }else{
                    if(weighted == 0){
                        writeLine("(1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} (max((EUD0["+(o1.index +1)+"] - (sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j])),0))^2) \n", bwParametersFile); /*Maicholl*/
                        writeLine("+ (1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} (max((sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j]) - OAR_targetUB,0))^2) \n", bwParametersFile); /*Maicholl*/

                    }else{
                        writeLine(""+wPar[count]+"*(1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} (max((EUD0["+(o1.index +1)+"] - (sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j])),0))^2) \n", bwParametersFile); /*Maicholl*/
                        writeLine("+ "+wPar[count]+"*(1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} (max((sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j]) - OAR_targetUB,0))^2) \n", bwParametersFile); /*Maicholl*/

                    }
                }
                count++;
            }
           /* for (Organs o1 : o) {
                if (!o1.isTarget){
                   //  writeLine("+ w["+(o1.index+1)+"] * ((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} (sum {j in 1..bmlt} "
                            + "x[j]*ddm"+o1.name+"[i,j])^a["+(o1.index +1)+"]))^(1/a["+(o1.index +1)+"])", bwParametersFile);
                   // 
                   writeLine("+ (1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} (max((sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j]) - EUD0["+(o1.index +1)+"],0))^4) \n", bwParametersFile); //Maicholl/
                }else{
                    writeLine("(1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} (max((EUD0["+(o1.index +1)+"] - (sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j])),0))^4) \n", bwParametersFile); /Maicholl/
                    //writeLine("+ (1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} (max((sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j]) - EUD0["+(o.length +1)+"],0)^2)) \n", bwParametersFile); /Maicholl/
                }
            }
            */
            writeLine(";\n s.t. \n", bwParametersFile);

            for (Organs o1 : o) {
                if (o1.isTarget){
                    writeLine("#equalityTarget: 	((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} (sum {j in 1..bmlt} "
                            + "x[j]*ddm"+o1.name+"[i,j])^a["+(o1.index +1)+"]))^(1/a["+(o1.index +1)+"]) >= t; \n", bwParametersFile);
                    writeLine("#constraintOAR_Target: 	((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} (sum {j in 1..bmlt} "
                            + "x[j]*ddm"+o1.name+"[i,j])^a["+(o.length +1)+"]))^(1/a["+(o.length +1)+"]) <=OAR_targetUB;\n", bwParametersFile);
                }else{
                    writeLine("#OAR_UB"+(o1.index +1)+": 	((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} (sum {j in 1..bmlt} "
                            + "x[j]*ddm"+o1.name+"[i,j])^a["+(o1.index +1)+"]))^(1/a["+(o1.index +1)+"]) <= UB"+(o1.index +1)+"; \n", bwParametersFile);
                    writeLine("MaxDosisM"+(o1.index +1)+": 	((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} (sum {j in 1..bmlt} "
                            + "x[j]*ddm"+o1.name+"[i,j]))) <= (EUD0["+(o1.index +1)+"]/2); \n", bwParametersFile);
                }
            }
            bwParametersFile.close();
        }
        
        
        System.out.println("Generating "+this.jobThreadID + "scriptCuadraticSum.sh");
        //Running Process
        scriptFile = this.jobThreadID + "scriptCuadraticSum.sh";
        System.out.println("Calling AMPL("+this.jobThreadID +")");
        Process p = new ProcessBuilder("ampl", scriptFile).start();
        InputStream is = p.getInputStream();
        InputStreamReader isr = new InputStreamReader(is);
        BufferedReader br = new BufferedReader(isr);
        String line;
        System.out.print("Solving CuadraticSum Function: ");
        //System.out.print("Angles: ");
        for(int i=0; i< this.beams; i++){
            System.out.print(" " +this.angles[i]+ " -- ");
        }
        System.out.print(" // ");
        for(int i=0; i< this.EUD0Par.length; i++){
            System.out.print(this.EUD0Par[i]+ " - ");
        }
        System.out.println();
        
        while ((line = br.readLine()) != null) {
            System.out.println(line);
        }
        is.close();
        isr.close();
        br.close();
        System.out.println("AMPL Process ("+this.jobThreadID +") has been destroyed");
        p.destroy();
        
    }
        
    
    
    public void deleteSolutionFiles(Organs[] o) throws IOException{
        System.out.println("Deleting Solutions Files from previous iterations for JobID "+jobThreadID+" ....Done");
        try{
            File file = new File("./"+this.jobThreadID + "currentSol.txt");
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(5)");
    		}
            }
            for (Organs o1 : o) {
                file = new File("./"+this.jobThreadID + "DVH_" + o1.name + ".txt");
                if (file.exists()) {
                    if(!file.delete()){
                        System.out.println("Delete operation for ./"+this.jobThreadID + "DVH_" + o1.name + ".txt failed.");
                    }
                }
                file = new File("./"+this.jobThreadID + "DVH_" + o1.name + ".txt");
                if (file.exists()) {
                    if(!file.delete()){
                        System.out.println("Delete operation for ./"+this.jobThreadID + "gEUD" + o1.name + ".txt failed.");
                    }
                }    
                if (o1.isTarget){
                    if (file.exists()) {
                        if(!file.delete()){
                            System.out.println("Delete operation for ./"+this.jobThreadID + "gEUD" + o1.name + "_UB.txt failed.");
                        }
                    }
                }
                if (file.exists()) {
                    if(!file.delete()){
                        System.out.println("Delete operation for ./"+this.jobThreadID + "voxelIndex" + o1.name + ".txt failed.");
                    }
                }
            }
    	}catch(Exception e){
    	}
        
    }
    
    public void generateParametersFile_cuadraticSum() throws IOException{ //Maicholl
        String parameterFile = "./"+this.jobThreadID + "extraCuadraticSum.dat";
        //Deleting existing parameter file extraCuadraticSum.dat
        File file = new File(parameterFile);
        if (file.exists()) {
            if(!file.delete()){
                System.out.println("Delete operation failed.(10)");
            }
        }
        //creating the new file
        
        BufferedWriter bwParametersFile=null;
        bwParametersFile =createBufferedWriter(parameterFile);
        
        writeLine("param a := ", bwParametersFile);
        for (int i=0;i<this.organs;i++){
            int j = i+1;
            writeLine(j + " " + this.aPar[i] + "\t", bwParametersFile);
        }
        //OAR-Target
        writeLine((this.organs+1) + " 10\t", bwParametersFile);
        writeLine(";\n", bwParametersFile);
        
        //writeLine(";\n", bwParametersFile);
        
        writeLine("#param w := ", bwParametersFile);
        for (int i=0;i<this.organs;i++){
            int j = i+1;
            writeLine(j + " " + this.wPar[i] + "\t", bwParametersFile);
        }
        //OAR-Target
        writeLine((this.organs+1) +" 8\t", bwParametersFile);
        writeLine(";\n", bwParametersFile);
        
        writeLine("param EUD0 := ", bwParametersFile);
        for (int i=0;i<this.organs;i++){
            int j = i+1;
            writeLine(j + " " + this.EUD0Par[i] + "\t", bwParametersFile);
        }
        writeLine(";\n", bwParametersFile);
        
        for (int i=0;i<this.organs;i++){
            int j=i+1;
            writeLine("param R" + j + " := " + this.voxels[i] + ";\n", bwParametersFile);
            
            if (this.isTarget[i]){
                writeLine("param t := " + EUD0Par[i] + ";\n", bwParametersFile);
                writeLine("param OAR_targetUB := " + this.UB[i] + ";\n", bwParametersFile);
            }else{
                writeLine("param UB" + j + " := " + this.UB[i] + ";\n", bwParametersFile);
                //writeLine("param LB" + j + " := " + this.LB[i] + ";\n", bwParametersFile);
            }
        }
        writeLine("param bmlt := " + this.totalBmlts + ";\n", bwParametersFile);
        writeLine("#param epsilon := " + this.epsilon + ";\n", bwParametersFile);
        
        bwParametersFile.close();  
        
    }
    
    public void generateParametersFile_cuadraticSum(double[] x) throws IOException{ //Maicholl
        String parameterFile = "./"+this.jobThreadID + "extraCuadraticSum.dat";
        //Deleting parameter file extra.txt
        
        System.out.println(parameterFile);  
        try{
            File file = new File(parameterFile);
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(3)");
    		}
            }
    	}catch(Exception e){
            e.printStackTrace();
    	}
        Random r = new Random();
        //creating the new file
        
        BufferedWriter bwParametersFile=null;
        bwParametersFile =createBufferedWriter(parameterFile);
        
        writeLine("var x := ", bwParametersFile);
        for (int i=0;i<this.totalBmlts;i++){
            int j = i+1;
            writeLine(j + " " + x[i] + "\t", bwParametersFile);
        }
        writeLine(";\n", bwParametersFile);
        
        writeLine("param a := ", bwParametersFile);
        for (int i=0;i<this.organs;i++){
            int j = i+1;
            writeLine(j + " " + this.aPar[i] + "\t", bwParametersFile);
        }
        //OAR-Target
        writeLine((this.organs+1) + " 60\t", bwParametersFile);
        writeLine(";\n", bwParametersFile);
        
        //writeLine(";\n", bwParametersFile);
        
        writeLine("#param w := ", bwParametersFile);
        for (int i=0;i<this.organs;i++){
            int j = i+1;
            writeLine(j + " " + this.wPar[i] + "\t", bwParametersFile);
        }
        //OAR-Target
        writeLine((this.organs+1) +" 8\t", bwParametersFile);
        writeLine(";\n", bwParametersFile);
        
        writeLine("param EUD0 := ", bwParametersFile);
        for (int i=0;i<this.organs;i++){
            int j = i+1;
            writeLine(j + " " + this.EUD0Par[i] + "\t", bwParametersFile);
            /*if(i+1==this.organs){
                writeLine(j+1 + " " + 85.8 + "\t", bwParametersFile);//Necesita el EUD del PTV como OAR Maicholl
            }     */  
        }
        writeLine(";\n", bwParametersFile);
        
        for (int i=0;i<this.organs;i++){
            int j=i+1;
            writeLine("param R" + j + " := " + this.voxels[i] + ";\n", bwParametersFile);
            
            if (this.isTarget[i]){
                writeLine("param t := " + EUD0Par[i] + ";\n", bwParametersFile);
                writeLine("param OAR_targetUB := " + this.UB[i] + ";\n", bwParametersFile);
            }else{
                writeLine("param UB" + j + " := " + this.UB[i] + ";\n", bwParametersFile);
                //writeLine("param LB" + j + " := " + this.LB[i] + ";\n", bwParametersFile);
            }
        }
        writeLine("param bmlt := " + this.totalBmlts + ";\n", bwParametersFile);
        writeLine("#param epsilon := " + this.epsilon + ";\n", bwParametersFile);
        
        bwParametersFile.close();
        
    }
    
    @SuppressWarnings("ConvertToStringSwitch")
    public void runLexBladder_Solver(Organs[] o) throws IOException{
        //Deleting Solution File currentSol.txt
        try{
            File file = new File(this.jobThreadID + "currentSol.txt");
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(10)");
    		}
            }
    	}catch(Exception e){
    	}
        
        //Create scriptLexicoRectum.sh file
        System.out.println("Creating " + this.jobThreadID + "scriptLexicoBladder.sh");  
        String scriptFile = this.jobThreadID + "scriptLexicoBladder.sh";
        try{
            File file = new File(scriptFile);
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(10.1)");
    		}
            }
    	}catch(Exception e){
    	}
        BufferedWriter bwParametersFile=null;
        bwParametersFile =createBufferedWriter(scriptFile);
        writeLine("model "+this.jobThreadID+"lexicoBladderModel.mod;\n", bwParametersFile);
        for (Organs o1 : o) {
            writeLine("data "+ this.jobThreadID + "DDM_" + o1.name + ".dat;\n", bwParametersFile);
        }
        
        writeLine("data "+ this.jobThreadID +"extra.dat;\n", bwParametersFile);
        writeLine("solve;\n", bwParametersFile);
        writeLine("display x > "+ this.jobThreadID +"currentSol.txt;\n", bwParametersFile);

        for (Organs o1 : o) {
            writeLine("display ((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} "
                        + "(sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j])^a["+(o1.index +1)+"]))^(1/a["+(o1.index +1)+"])" + " > " + this.jobThreadID + "gEUD_" + o1.name + ".txt;\n", bwParametersFile); 
            writeLine(" display {i in 1 .. R"+(o1.index +1)+"} (sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j]) > " + this.jobThreadID + "DVH_" + o1.name + ".txt;\n", bwParametersFile); 
            if (o1.isTarget){
                writeLine("display " + "((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} "
                        + "(sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j])^a["+(o.length +1)+"]))^(1/a["+(o.length +1)+"]) > " + this.jobThreadID + "gEUD_" + o1.name + "_UB.txt;\n", bwParametersFile); 
            }
        }
        writeLine("display _con > "+ this.jobThreadID +"lagrangeMultiplier.txt;\n", bwParametersFile);
        bwParametersFile.close();
        
        //Create LexicoBladder.mod for AMPL*/
        System.out.println("Creating " + this.jobThreadID + "lexicoBladderModel.mod");  
        scriptFile = this.jobThreadID + "lexicoBladderModel.mod";
        try{
            File file = new File(scriptFile);
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(10.2)");
    		}
            }
    	}catch(Exception e){
    	}
        //creating the new file
        bwParametersFile=null;
        bwParametersFile =createBufferedWriter(scriptFile);
        
        if (null != solver)switch (solver) {
            case "ipopt":
                writeLine("option solver ipopt; \n", bwParametersFile);
                //writeLine("options ipopt_options \"linear_solver=ma57 linear_system_scaling=mc19 wantsol=8 print_level=4 tol=0.0001\"; \n", bwParametersFile);
                writeLine("options ipopt_options \"wantsol=8 print_level=0 tol=0.0001\"; \n", bwParametersFile);
                break;
            case "knitro":
                System.out.println("Entro knitro. Solver: " + solver +". Function: 13");
                writeLine("option solver knitro; \n", bwParametersFile);
                writeLine("option knitro_options \"wantsol=8 outlev=1\"; \n", bwParametersFile);
                break;
            default:
                writeLine("option solver ipopt; \n", bwParametersFile);
                //writeLine("options ipopt_options \"linear_solver=ma57 linear_system_scaling=mc19 wantsol=8 print_level=4 tol=0.0001\"; \n", bwParametersFile);
                writeLine("options ipopt_options \"wantsol=8 print_level=0 tol=0.0001\"; \n", bwParametersFile);
                
        }
        
        for (Organs o1 : o) {
            writeLine("param R" + (o1.index +1) + "; #number of voxels of " + o1.name + "\n", bwParametersFile);
        }
        writeLine("param bmlt; 	#number of beamlets \n", bwParametersFile);
        
        for (Organs o1 : o) {
            writeLine("param ddm" + o1.name + "{1 .. R" + (o1.index +1) + ", 1 .. bmlt};\n", bwParametersFile);
        }
        for (Organs o1 : o) {
            writeLine("param UB" + (o1.index +1) + ";\n", bwParametersFile);
        }
        writeLine("param t;\n" + "param epsilon;\n" + "param OAR_targetUB;\n", bwParametersFile);
        writeLine("var x {1 .. bmlt} >= 0, <=20, default 1; \n", bwParametersFile);
        
        writeLine("param a{1 .. "+ (o.length + 1) + "}; \n", bwParametersFile);
        
        writeLine("minimize Total_Cost: ", bwParametersFile);
        
        for (Organs o1 : o) {
            if ("BLADDER".equals(o1.name)|| "VEJIGA".equals(o1.name)){
                writeLine("((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} "
                        + "(sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j])^a["+(o1.index +1)+"]))^(1/a["+(o1.index +1)+"]);\n", bwParametersFile); 
            }
        }
        writeLine(";\n s.t. \n", bwParametersFile);
        
        for (Organs o1 : o) {
            if (o1.isTarget){
                writeLine("equalityTarget: 	((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} (sum {j in 1..bmlt} "
                        + "x[j]*ddm"+o1.name+"[i,j])^a["+(o1.index +1)+"]))^(1/a["+(o1.index +1)+"]) >= t; \n", bwParametersFile);
                writeLine("constraintOAR_Target: 	((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} (sum {j in 1..bmlt} "
                        + "x[j]*ddm"+o1.name+"[i,j])^a["+(o.length +1)+"]))^(1/a["+(o.length +1)+"]) <=OAR_targetUB;\n", bwParametersFile);
            }
        }
        bwParametersFile.close();
            
        
        //Running Process
        scriptFile = this.jobThreadID + "scriptLexicoBladder.sh";
        Process p = new ProcessBuilder("ampl", scriptFile).start();
        InputStream is = p.getInputStream();
        InputStreamReader isr = new InputStreamReader(is);
        BufferedReader br = new BufferedReader(isr);
        String line;
        System.out.print("Solving LexicoBladder problem: ");
        //System.out.println("Angles: ");
        for(int i=0; i< this.beams; i++){
            System.out.print(" " +this.angles[i]+ " -- ");
        }
        System.out.print(" // ");
        for(int i=0; i< this.EUD0Par.length; i++){
            System.out.print(this.EUD0Par[i]+ " - ");
        }
        System.out.println();
        while ((line = br.readLine()) != null) {
            System.out.println(line);
        }
        isr.close();
        is.close();
        br.close();
        p.destroy();
    }
    
    public void writeLine(String l, BufferedWriter bw) throws IOException{
	
	String row;
	row = l; 
	bw.write(row);
	

    }
   
}
