################################################################################
# CQRNDs2_equalsize: Circular Quasi Rees neighbor designs for block of equal size(K)
################################################################################
# Algorithm from paper:
# Akbar Firdos1, Mahmood Ul Hassan2, Farrukh Jamal1, Hurria Ali1, Khadija Noreen1 and Rashid Ahmed1. Efficient Balanced and Strongly 
# Balanced Neighbor Designs through  CQRNDs2 -designs Generated with R 
# Coded by Ali et al., 01-03-2022 to 05-07-2022
# Version 2.0  (2022-07-05)
################################################################
# Division of adjusted A in i-1 groups of size k and one
# of size K-2 to get the set(s) of shifts
################################################################


Designs<-function(v,k1,k2,k3,i){
  if(v>0 & k1>0 &k2==0 &k3==0){
    ################################################################################
    # CQRNDs2_equalsize: Circular Quasi Rees neighbor designs for block of equal size(K)
    ################################################################################
    # Algorithm from paper:
    # Akbar Firdos1, Mahmood Ul Hassan2, Farrukh Jamal1, Hurria Ali1, Khadija Noreen1 and Rashid Ahmed1. Efficient Balanced and Strongly
    # Balanced Neighbor Designs through  CQRNDs2 -designs Generated with R
    # Coded by Ali et al., 01-03-2022 to 05-07-2022
    # Version 2.0  (2022-07-05)
    ################################################################
    # Division of adjusted A in i-1 groups of size k and one
    # of size K-2 to get the set(s) of shifts
    ################################################################
    grouping1<-function(A,k,v,i){
      bs<-c()
      z=0;f=1
      A1=A
      while(f<=((i+1)-1)){
        for(y in 1:5000){
          comp<-sample(1:length(A1),k)
          com<-A1[comp]
          cs<-sum(com)
          if(cs%%v==0){
            bs<-rbind(bs,com)
            A1<-A1[-comp]
            z<-z+1
            f=f+1
          }
          if(z==i) break
        }
        if(z<i) {bs<-c();z=0;f=1;A1=A}
      }
      bs1<-t(apply(bs,1,sort))
      bs1<-cbind(bs1,rowSums(bs),rowSums(bs)/v)
      rownames(bs1)<-paste("G",1:i, sep="")
      colnames(bs1)<-c(paste(1:k, sep=""),"sum" ,"sum/v")
      bs2<-t(apply(bs,1,sort))
      bs2<-delmin(bs2)
      list(B1=list(bs2),B2=list(bs1),B3=A1)
    }
    #######################################################################
    # Obtaing set(s) of shifts by deleting smallest value of each group
    #######################################################################
    delmin<-function(z){
      fs<-c()
      n<-nrow(z)
      c<-ncol(z)-1
      for(i in 1:n){
        z1<-z[i,]
        z2<-z1[z1!=min(z1)]
        fs<-rbind(fs,z2)
      }
      rownames(fs)<-paste("S",1:n, sep="")
      colnames(fs)<-rep("",c)
      return(fs)
    }
    ################################################################################
    # Selection of adjusted A and the set(s) of shifs to obtain Circular
    # balance neighbor design for block of equal size.
    ################################################################################
    # D=1: Circular Balanced Neighbor Designs
    # D=2: Circular Strongly Balanced Neighbor Designs
    #   K: Block sizes
    #   i: Number of set of shifts for K
    CBND_equalsize<-function(k,i,D=1){
      if(k<=2) stop("k= Block size: Block size must be greater than 2")
      setClass( "stat_test", representation("list"))
      setMethod("show", "stat_test", function(object) {
        row <- paste(rep("", 70), collapse = "")
        cat(row, "\n")
        row <- paste(rep("=", 40), collapse = "")
        cat(row, "\n")
        print(object$S[[1]])
      })
      if(D==2){
        v=2*i*k-1; m=(v-1)/2
        if(m%%8==0){
          j=m/8
          if(j<1) {return("Conditions are not satisfied for CSBND")}
          A=c(0:(j-1),(j+1):m,(v-j))
          A1<-grouping1(A,k,v,i)
          A2<-c(v,k);names(A2)<-c("V","K")
          x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
        }
        if(m%%8==1){
          j=(m-1)/8
          if(j<1) {return("Conditions are not satisfied for CSBND")}
          A=c(0:(3*j),(3*j+2):(m-1),(m+1),(v-(3*j+1)))
          A1<-grouping1(A,k,v,i)
          A2<-c(v,k);names(A2)<-c("V","K")
          x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
        }
        if(m%%8==2){
          j=(m-2)/8
          if(j<1) {return("Conditions are not satisfied for CSBND")}
          A=c(0:(5*j+1),(5*j+3):(m-1),(m+1),(v-(5*j+2)))
          A1<-grouping1(A,k,v,i)
          A2<-c(v,k);names(A2)<-c("V","K")
          x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
        }
        if(m%%8==3){
          j=(m-3)/8
          if(j<0) {return("Conditions are not satisfied for CSBND")}
          A=c(0:(m-j-1),(m-j+1):m,(v-(m-j)))
          A1<-grouping1(A,k,v,i)
          A2<-c(v,k);names(A2)<-c("V","K")
          x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
        }
        if(m%%8==4){
          j=(m-4)/8
          if(j<0) {return("Conditions are not satisfied for CSBND")}
          A=c(0:j,(j+2):(m-1),(m+1),(v-(j+1)))
          A1<-grouping1(A,k,v,i)
          A2<-c(v,k);names(A2)<-c("V","K")
          x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
        }
        if(m%%8==5){
          j=(m-5)/8
          if(j<0) {return("Conditions are not satisfied for CSBND")}
          A=c(0:(3*j+1),(3*j+3):(m),(v-(3*j+2)))
          A1<-grouping1(A,k,v,i)
          A2<-c(v,k);names(A2)<-c("V","K")
          x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
        }
        if(m%%8==6){
          j=(m-6)/8
          if(j<0) {return("Conditions are not satisfied for CSBND")}
          A=c(0:(5*j+3),(5*j+5):(m),(v-(5*j+4)))
          A1<-grouping1(A,k,v,i)
          A2<-c(v,k);names(A2)<-c("V","K")
          x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
        }
        if(m%%8==7){
          j=(m-7)/8
          if(j<1) {return("Conditions are not satisfied for CSBND")}
          A=c(0:(m-j-1),(m-j+1):(m-1),(m+1),(v-(m-j)))
          A1<-grouping1(A,k,v,i)
          A2<-c(v,k);names(A2)<-c("V","K")
          x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
        }  }
      if(D==1){
        v=2*i*k+1; m=(v-1)/2
        if(m%%8==0){
          j=m/8
          if(j<1) {return("Conditions are not satisfied for CBNDs")}
          A=c(1:(j-1),(j+1):m,(v-j))
          A1<-grouping1(A,k,v,i)
          A2<-c(v,k);names(A2)<-c("V","K")
          x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
        }
        if(m%%8==1){
          j=(m-1)/8
          if(j<1) {return("Conditions are not satisfied for CBNDs")}
          A=c(1:(3*j),(3*j+2):(m-1),(m+1),(v-(3*j+1)))
          A1<-grouping1(A,k,v,i)
          A2<-c(v,k);names(A2)<-c("V","K")
          x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
        }
        if(m%%8==2){
          j=(m-2)/8
          if(j<1) {return("Conditions are not satisfied for CBNDs")}
          A=c(1:(5*j+1),(5*j+3):(m-1),(m+1),(v-(5*j+2)))
          A1<-grouping1(A,k,v,i)
          A2<-c(v,k);names(A2)<-c("V","K")
          x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
        }
        if(m%%8==3){
          j=(m-3)/8
          if(j<0) {return("Conditions are not satisfied for CBNDs")}
          A=c(1:(m-j-1),(m-j+1):m,(v-(m-j)))
          A1<-grouping1(A,k,v,i)
          A2<-c(v,k);names(A2)<-c("V","K")
          x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
        }
        if(m%%8==4){
          j=(m-4)/8
          if(j<0) {return("Conditions are not satisfied for CBNDs")}
          A=c(1:j,(j+2):(m-1),(m+1),(v-(j+1)))
          A1<-grouping1(A,k,v,i)
          A2<-c(v,k);names(A2)<-c("V","K")
          x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
        }
        if(m%%8==5){
          j=(m-5)/8
          if(j<0) {return("Conditions are not satisfied for CBNDs")}
          A=c(1:(3*j+1),(3*j+3):(m),(v-(3*j+2)))
          A1<-grouping1(A,k,v,i)
          A2<-c(v,k);names(A2)<-c("V","K")
          x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
        }
        if(m%%8==6){
          j=(m-6)/8
          if(j<0) {return("Conditions are not satisfied for CBNDs")}
          A=c(1:(5*j+3),(5*j+5):(m),(v-(5*j+4)))
          A1<-grouping1(A,k,v,i)
          A2<-c(v,k);names(A2)<-c("V","K")
          x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
        }
        if(m%%8==7){
          j=(m-7)/8
          if(j<1) {return("Conditions are not satisfied for CBNDs")}
          A=c(1:(m-j-1),(m-j+1):(m-1),(m+1),(v-(m-j)))
          A1<-grouping1(A,k,v,i)
          A2<-c(v,k);names(A2)<-c("V","K")
          x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
        }
      }
      if(D==12){
        v=2*i*k-3; m=(v-1)/2
        if(m%%8==0){
          j=m/8
          if(j<1) {return("Conditions are not satisfied for CSBGND-Odd")}
         
           A=c(0:(j-1),(j+1),(j+2):(m-1),(m+1),(2*m-1-j),1)
          A1<-grouping1(A,k,v,i)
          A2<-c(v,k);names(A2)<-c("V","K")
          x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
        }
        if(m%%8==1){
          j=(m-1)/8
          if(j<1) {return("Conditions are not satisfied for CSBGND-Odd")}
          A=c(0:(3*j),(3*j+2),(3*j+3):m,(2*m-3*j),1)
          A1<-grouping1(A,k,v,i)
          A2<-c(v,k);names(A2)<-c("V","K")
          x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
        }
        if(m%%8==2){
          j=(m-2)/8
          if(j<1) {return("Conditions are not satisfied for CSBGND-Odd")}
          A=c(0:(5*j+1),(5*j+3),(5*j+4):m,(2*m-5*j-1),1)
          A1<-grouping1(A,k,v,i)
          A2<-c(v,k);names(A2)<-c("V","K")
          x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
        }
        if(m%%8==3){
          j=(m-3)/8
          if(j<0) {return("Conditions are not satisfied for CSBGND-Odd")}
          A=c(0:(m-j),(m-j+2),(m-j+3):m,(m+j),1)
          A1<-grouping1(A,k,v,i)
          A2<-c(v,k);names(A2)<-c("V","K")
          x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
        }
        if(m%%8==4){
          j=(m-4)/8
          if(j<0) {return("Conditions are not satisfied for CSBGND-Odd")}
          A=c(0:j,(j+2),(j+3):m,(2*m-j),1)
          A1<-grouping1(A,k,v,i)
          A2<-c(v,k);names(A2)<-c("V","K")
          x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
        }
        if(m%%8==5){
          j=(m-5)/8
          if(j<0) {return("Conditions are not satisfied for CSBGND-Odd")}
          A=c(0:(3*j+2),(3*j+4),(3*j+5):(m-1),(m+1),(2*m-3*j-2),1)
          A1<-grouping1(A,k,v,i)
          A2<-c(v,k);names(A2)<-c("V","K")
          x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
        }
        if(m%%8==6){
          j=(m-6)/8
          if(j<0) {return("Conditions are not satisfied for CSBGND-Odd")}
          A=c(0:(5*j+4),(5*j+6),(5*j+7):(m-1),(m+1),(2*m-5*j-4),1)
          A1<-grouping1(A,k,v,i)
          A2<-c(v,k);names(A2)<-c("V","K")
          x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
        }
        if(m%%8==7){
          j=(m-7)/8
          if(j<1) {return("Conditions are not satisfied for CSBGND-Odd")}
          A=c(0:(m-j-1),(m-j+1),(m+j+2):m,(m+j+1),1)
          A1<-grouping1(A,k,v,i)
          A2<-c(v,k);names(A2)<-c("V","K")
          x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
        }    }
      if(D==11){
        v=2*i*k+3; m=(v-1)/2
        if(m%%8==0){
          j=m/8
          if(j<1) {return("Conditions are not satisfied for CPBND-Odd")}
          A=c(2:j,(j+2),(j+3):(m-2),m,(m+2),(2*m-j))
          A1<-grouping1(A,k,v,i)
          A2<-c(v,k);names(A2)<-c("V","K")
          x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
        }
        if(m%%8==1){
          j=(m-1)/8
          if(j<1) {return("Conditions are not satisfied for CPBND-Odd")}
          A=c(2:(3*j-1),(3*j+1),(3*j+2):m,(2*m-(3*j-1)))
          A1<-grouping1(A,k,v,i)
          A2<-c(v,k);names(A2)<-c("V","K")
          x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
        }
        if(m%%8==2){
          j=(m-2)/8
          if(j<1) {return("Conditions are not satisfied for CPBND-Odd")}
          A=c(2:(5*j+2),(5*j+3):m,(2*m-(5*j)))
          A1<-grouping1(A,k,v,i)
          A2<-c(v,k);names(A2)<-c("V","K")
          x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
        }
        if(m%%8==3){
          j=(m-3)/8
          if(j<0) {return("Conditions are not satisfied for CPBND-Odd")}
          A=c(2:(m-j-1),(m-j+1),(m-j+2):(m+1),(m+2),(m+j+1))
          A1<-grouping1(A,k,v,i)
          A2<-c(v,k);names(A2)<-c("V","K")
          x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
        }
        if(m%%8==4){
          j=(m-4)/8
          if(j<0) {return("Conditions are not satisfied for CPBND-Odd")}
          A=c(2:(j+1),(j+3),(j+4):(m-2),(m+1),(m+2),(2*m-j-1))
          A1<-grouping1(A,k,v,i)
          A2<-c(v,k);names(A2)<-c("V","K")
          x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
        }
        if(m%%8==5){
          j=(m-5)/8
          if(j<0) {return("Conditions are not satisfied for CPBND-Odd")}
          A=c(2:(3*j+1),(3*j+3),(3*j+4):(m-1),(m+1),(2*m-3*j-1))
          A1<-grouping1(A,k,v,i)
          A2<-c(v,k);names(A2)<-c("V","K")
          x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
        }
        if(m%%8==6){
          j=(m-6)/8
          if(j<0) {return("Conditions are not satisfied for CPBND-Odd")}
          A=c(2:(5*j+3),(5*j+5),(5*j+6):(m-1),(m+1),(2*m-5*j-3))
          A1<-grouping1(A,k,v,i)
          A2<-c(v,k);names(A2)<-c("V","K")
          x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
        }
        if(m%%8==7){
          j=(m-7)/8
          if(j<1) {return("Conditions are not satisfied for CPBND-Odd")}
          A=c(2:(m-j-2),(m-j),(m-j+1):m,(m+j+2))
          A1<-grouping1(A,k,v,i)
          A2<-c(v,k);names(A2)<-c("V","K")
          x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
        }  }
      if(D==3){
        v=2*i*k+2; m=(v-2)/2
        if(m%%4==0){
          A=1:m
          A1<-grouping1(A,k,v,i)
          A2<-c(v,k);names(A2)<-c("V","K")
          x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
        }
        if(m%%4==3){
          A<-c(1:((3*m-1)/4),((3*m+7)/4),((3*m+11)/4):m, (5*(m+1)/4))
          A1<-grouping1(A,k,v,i)
          A2<-c(v,k);names(A2)<-c("V","K")
          x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
        }
        if(m%%4==1 |  m%%4==2){return("The minimal CPNDs in which v/2 unordered pairs cannot be constructed for v=2ik+2 and k=block size")}
      }
      if(D==4){
        v=2*i*k+4; m=(v-2)/2
        if(m%%4==0){
          A=c(1:(m/2),((m+4)/2),((m+6)/2):(m-1),((3*m+2)/2))
          A1<-grouping1(A,k,v,i)
          A2<-c(v,k);names(A2)<-c("V","K")
          x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
        }
        if(m%%4==1){
          A=c(1:((3*m+1)/4),((3*m+9)/4),((3*m+13)/4):(m-1),((5*m+3)/4))
          A1<-grouping1(A,k,v,i)
          A2<-c(v,k);names(A2)<-c("V","K")
          x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
        }
        if(m%%4==2){
          A=c(2:(m-2),m,(2*m+1))
          A1<-grouping1(A,k,v,i)
          A2<-c(v,k);names(A2)<-c("V","K")
          x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
        }
        if(m%%4==3){
          A=c(1:((m+1)/4),((m+9)/4),((m+13)/4):(m-2), m, ((7*m+3)/4))
          A1<-grouping1(A,k,v,i)
          A2<-c(v,k);names(A2)<-c("V","K")
          x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
        }  }
      if(D==5){
        v=2*i*k; m=(v-2)/2
        if(m%%4==2){
          A=1:(m+1)
          A1<-grouping1(A,k,v,i)
          A2<-c(v,k);names(A2)<-c("V","K")
          x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
        }
        if(m%%4==3){
          A<-c(1:((m-3)/4),((m+5)/4),((m+9)/4):m,(m+1),(7*(m+1)/4))
          A1<-grouping1(A,k,v,i)
          A2<-c(v,k);names(A2)<-c("V","K")
          x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
        }
        if(m%%4==1 |  m%%4==0){return("The minimal CWBNDs in which v/2 unordered pairs apears twice cannot be constructed for v=2ik and k=block size")}
      }
      if(D==6){
        v=2*i*k-2; m=(v-2)/2
        if(m%%4==0){
          A=c(1:(m-1),(m-1),(m+1),(m+2))
          A1<-grouping1(A,k,v,i)
          A2<-c(v,k);names(A2)<-c("V","K")
          x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
        }
        if(m%%4==1){
          A=c(1:((m-5)/4),((m+3)/4),((m+7)/4):m,(m+1),m,((7*m+9)/4))
          A1<-grouping1(A,k,v,i)
          A2<-c(v,k);names(A2)<-c("V","K")
          x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
        }
        if(m%%4==2){
          A=c(1:((m-2)/2),((m+2)/2),((m+4)/2):m,(m+1),m,((3*m+4)/2))
          A1<-grouping1(A,k,v,i)
          A2<-c(v,k);names(A2)<-c("V","K")
          x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
        }
        if(m%%4==3){
          A=c(1:((3*m-5)/4),((3*m+3)/4),((3*m+7)/4):(m-1), m,(m+1),(m-1), ((5*m+9)/4))
          A1<-grouping1(A,k,v,i)
          A2<-c(v,k);names(A2)<-c("V","K")
          x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
        }  }
      if(D==7){
        v=2*i*k; m=(v-2)/2
        if(m%%4==0){
          A=0:m
          A1<-grouping1(A,k,v,i)
          A2<-c(v,k);names(A2)<-c("V","K")
          x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
        }
        if(m%%4==3){
          A<-c(0:((3*m-1)/4),((3*m+7)/4),((3*m+11)/4):m, (5*(m+1)/4))
          A1<-grouping1(A,k,v,i)
          A2<-c(v,k);names(A2)<-c("V","K")
          x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
        }
        if(m%%4==1 | m%%4==2){return("The minimal CSPBNDs in which v/2 unordered pairs cannot be constructed for v=2ik+2 and k=block size")}
      }
      if(D==8){
        v=2*i*k+2; m=(v-2)/2
        if(m%%4==0){
          A=c(0:(m/2),((m+4)/2),((m+6)/2):(m-1),((3*m+2)/2))
          A1<-grouping1(A,k,v,i)
          A2<-c(v,k);names(A2)<-c("V","K")
          x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
        }
        if(m%%4==1){
          A=c(0:((3*m+1)/4),((3*m+9)/4),((3*m+13)/4):(m-1),((5*m+3)/4))
          A1<-grouping1(A,k,v,i)
          A2<-c(v,k);names(A2)<-c("V","K")
          x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
        }
        if(m%%4==2){
          A=c(0,2:(m-2),m,(2*m+1))
          A1<-grouping1(A,k,v,i)
          A2<-c(v,k);names(A2)<-c("V","K")
          x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
        }
        if(m%%4==3){
          A=c(0:((m+1)/4),((m+9)/4),((m+13)/4):(m-2), m, ((7*m+3)/4))
          A1<-grouping1(A,k,v,i)
          A2<-c(v,k);names(A2)<-c("V","K")
          x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
        }
      }
      if(D==9){
        v=(2*i*k)-2; m=(v-2)/2
        if(m%%4==2){
          A=0:(m+1)
          A1<-grouping1(A,k,v,i)
          A2<-c(v,k);names(A2)<-c("V","K")
          x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
        }
        if(m%%4==3){
          A<-c(0:((m-3)/4),((m+5)/4),((m+9)/4):m,(m+1),(7*(m+1)/4))
          A1<-grouping1(A,k,v,i)
          A2<-c(v,k);names(A2)<-c("V","K")
          x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
        }
        if(m%%4==1 |  m%%4==0){return("The minimal CSBGNDs in which v/2 unordered pairs apears twice cannot be constructed for v=2ik and k=block size")}
      }
      if(D==10){
        v=2*i*k-4; m=(v-2)/2
        if(m%%4==0){
          A=c(0:(m-1),(m-1),(m+1),(m+2))
          A1<-grouping1(A,k,v,i)
          A2<-c(v,k);names(A2)<-c("V","K")
          x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
        }
        if(m%%4==1){
          A=c(0:((m-5)/4),((m+3)/4),((m+7)/4):m,(m+1),m,((7*m+9)/4))
          A1<-grouping1(A,k,v,i)
          A2<-c(v,k);names(A2)<-c("V","K")
          x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
        }
        if(m%%4==2){
          A=c(0:((m-2)/2),((m+2)/2),((m+4)/2):m,(m+1),m,((3*m+4)/2))
          A1<-grouping1(A,k,v,i)
          A2<-c(v,k);names(A2)<-c("V","K")
          x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
        }
        if(m%%4==3){
          A=c(0:((3*m-5)/4),((3*m+3)/4),((3*m+7)/4):(m-1), m,(m+1),(m-1), ((5*m+9)/4))
          A1<-grouping1(A,k,v,i)
          A2<-c(v,k);names(A2)<-c("V","K")
          x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
        }  }
      new("stat_test", x)
    }
    design_CBND<-function(H){
      setClass( "CBND_design", representation("list"))
      setMethod("show", "CBND_design", function(object) {
        row <- paste(rep("=", 51), collapse = "")
        cat(row, "\n")
        cat("Following is minimal CBND for", "v=" ,object$R[1], "and","k=",object$R[2], "\n")
        row <- paste(rep("=", 51), collapse = "")
        cat(row, "\n")
        for(i in 1:length(ss)){
          W<-ss[[i]]
          nr<-dim(W)[1]
          for(j in 1:nr){
            print(object$Design[[i]][[j]])
            cat("\n\n")
          }}
      })
      v<-H$R[1]
      k<-H$R[2]
      ss<-H$S
      treat<-(1:v)-1
      fn<-(1:v)
      G<-list()
      for(j in 1:length(ss)){
        W<-ss[[j]]
        nr<-dim(W)[1]
        nc<-dim(W)[2]
        D<-list()
        for(i in 1:nr){
          dd<-c()
          d1<-matrix(treat,(nc+1),v,byrow = T)
          ss1<-cumsum(c(0,W[i,]))
          dd2<-d1+ss1
          dd<-rbind(dd,dd2)
          rr<-dd[which(dd>=v)]%%v
          dd[which(dd>=v)]<-rr
          colnames(dd)<-paste("B",fn, sep="")
          rownames(dd)<-rep("",(nc+1))
          fn<-fn+v
          D[[i]]<-dd
        }
        G[[j]]<-D
      }
      x<-list(Design=G,R=H$R)
      new("CBND_design", x)
    }
    Create_Designs<-function(v,k,i){
      if(v%%2!=0){
        if((v-1)%%(2*k)==0 & (i=(v-1)/(2*k)) &(i%%1==0)&(i>0)){
          cat("MCBNDs is possible for","v=",v,"k=",k,"i=",i,rep("", 30))
          print(H<-CBND_equalsize(k=k,i=i,D=1))
        }
          if((v+1)%%(2*k)==0 & (i=(v+1)/(2*k))&(i%%1==0)&(i>0) ) {
            cat("MCSBNDs is possible for","v=",v,"k=",k,"i=",i,rep("", 30))
            print(H<-CBND_equalsize(k=k,i=i,D=2))
          }
            if((v-3)%%(2*k)==0 & (i=(v-3)/(2*k))&(i%%1==0)&(i>0) ){
              cat("MCPBND-Odd is possible for","v=",v,"k=",k,"i=",i,rep("", 30))
              print(H<-CBND_equalsize(k=k,i=i,D=11))
            }
              if((v+3)%%(2*k)==0 & (i=(v+3)/(2*k)) &(i%%1==0) & (i>0)) {
                cat("MCSBGND-Odd is possible for","v=",v,"k=",k,"i=",i,rep("", 30))
                print(H<-CBND_equalsize(k=k,i=i,D=12))
              }
        if((i%%1<0) |(i%%1>0)) {cat("Other designs are not possible and change value ","v=",v,"and", "k=",k,"i=",i, " run again",rep("", 30))
}
          }else
        if(v%%2==0){
          if((v-2)%%(2*k)==0 & (i=(v-2)/(2*k))  & (((v-2)/2)%%4==0 |((v-2)/2)%%4==3)&(i%%1==3)&(i>0)){
            cat("MCPBNDs-I is possible for","v=",v,"k=",k,"i=",i,rep("",15))
            print(H<-CBND_equalsize(k=k,i=i,D=3))
          }
            if((v-2)%%(2*k)==0 & (i=(v-2)/(2*k))&(i%%1==2)&(i>0) ){
              cat(" MCSPBNDs-II is possible for","v=",v,"k=",k,"i=",i,rep("",15))
              print(H<-CBND_equalsize(k=k,i=i,D=8))
            }
              if((v-4)%%(2*k)==0 & (i=(v-4)/(2*k)) & (i%%1==0) & (i>0)){
                cat("MCPBNDs-II is possible for","v=",v,"k=",k,"i=",i,rep("",15))
                print(H<-CBND_equalsize(k=k,i=i,D=4))
              }
                
                  if(v%%(2*k)==0 & (i=v/(2*k)) &(i%%1==0)&(i>0)& (((v-2)/2)%%4==2 |((v-2)/2)%%4==3)){
                    cat("MCWBNDs-I  is possible for","v=",v,"k=",k,"i=",rep("",15))
                    print(H<-CBND_equalsize(k=k,i=i,D=5))
                  }
                    if(v%%(2*k)==0 & (i=v/(2*k)) &(i%%1==0)&(i>0)& (((v-2)/2)%%4==0 |((v-2)/2)%%4==3)){
                      cat("MCSPBNDs-I is possible for","v=",v,"k=",k,"i=",i,rep("",15))
                      print(H<-CBND_equalsize(k=k,i=i,D=7))
                    }
                      if((v+2)%%(2*k)==0 & (i=(v+2)/(2*k)) &(i%%1==0)&(i>0)& (((v-2)/2)%%4==2 |((v-2)/2)%%4==3)){
                        cat("MCSBGNDs-I is possible for","v=",v,"k=",k,"i=",i,rep("",15))
                        print(H<-CBND_equalsize(k=k,i=i,D=9))
                      }
                        if((v+2)%%(2*k)==0 & (i=(v+2)/(2*k))&(i%%1==0)&(i>0)){
                          cat("MCWBNDs-II  is possible for","v=",v,"k=",k,"i=",i,rep("",15))
                          print(H<-CBND_equalsize(k=k,i=i,D=6))
                        }
                          if((v+4)%%(2*k)==0 & (i=(v+4)/(2*k))&(i%%1==0)&(i>0)){
                            cat("MCSGNDs-II is possible for","v",v,"k=",k,"i=",i,rep("",15))
                            print(H<-CBND_equalsize(k=k,i=i,D=10))
                          }
          if((i%%1<0) |(i%%1>0)){cat("Other designs are not possible for given  value","v=",v,"and", "k=",k,"i=",i,rep("",15))
        
                          }  }  }
    ################################################################################
    # Examples: Using CBND_equal size function to obtain the set(s) of shifts
    # for construction of circular balance neighbor design for equal block
    # sizes (k)
    ################################################################################
    #Run the program to create designs and put the v and k value in the Create_Designs Function
    #Run this program get the possible designs
    (H<-Create_Designs(v=v, k=k1))
    
    
  }else
  if(v>0 &k1>0 &k2>0 &k3==0){################################################################################
    # CBND_2diffsize: Circular balance neighbor design for block of two different
    # sizes (K1 and k2)
    # Algorithm from paper:
    # Akbar Fardos, Khadija Noreen, Muhammad Sajid Rashid, Mahmood Ul Hassan, 
    # Zahra Noreen and Rashid Ahmed (2021). An Algorithm to Generate Minimal 
    # Circular Balanced and Strongly Balanced Neighbor Designs
    # Coded by Fardos et al., 01-08-2021 to 05-09-2021
    # Version 1.4.0  (2021-09-05)
    ################################################################################
    # Selection of i groups of size K1 from adjusted A. The set of remaining 
    # (Unselected) elements are saved in the object named as B2. 
    ################################################################################
    grouping1<-function(A,k,v,i){
      bs<-c()
      z=0;f=1
      A1=A
      while(f<=i){
        for(y in 1:5000){
          comp<-sample(1:length(A1),k)
          com<-A1[comp]
          cs<-sum(com)
          if(cs%%v==0){
            bs<-rbind(bs,com)
            A1<-A1[-comp]
            z<-z+1
            f=f+1
          }
          if(z==i) break
        }
        if(z<i) {bs<-c();z=0;f=1;A1=A}  
      }
      list(B1=bs,B2=A1)
    }
    ################################################################################
    # Selection of i group of size K1 from adjusted A and division of required 
    # number of groups of size K2 from B2. 
    ################################################################################
    grouping2<-function(A,k,v,i,sk2){
      bs1<-c()
      j=i+sk2
      z=0;f=1
      A1=A
      while(f<=j){
        s<-grouping1(A1,k[1],v,i)
        A2<-s$B2
        z=i;f=f+i
        for(y in 1:1000){
          comp<-sample(1:length(A2),k[2])
          com<-A2[comp]
          cs<-sum(com)
          if(cs%%v==0){
            bs1<-rbind(bs1,com)
            A2<-A2[-comp]
            z<-z+1
            f=f+1
          }
          if(z==j) break
        }
        if(z<j) {bs1<-c();z=0;f=1;A1=A}  
      }
      gs1<-t(apply(s$B1,1,sort))
      gs1<-cbind(gs1,rowSums(gs1),rowSums(gs1)/v)
      rownames(gs1)<-paste("G",1:i, sep="")
      colnames(gs1)<-c(paste(1:k[1], sep=""),"sum" ,"sum/v")
      gs2<-t(apply(bs1,1,sort))
      gs2<-cbind(gs2,rowSums(gs2),rowSums(gs2)/v)
      rownames(gs2)<-paste("G",(nrow(gs1)+1):(nrow(gs1)+sk2), sep="")
      colnames(gs2)<-c(paste(1:k[2], sep=""),"sum" ,"sum/v")
      fs1<-t(apply(s$B1,1,sort))
      fs1<-delmin(fs1)
      rownames(fs1)<-paste("S",1:i, sep="")
      colnames(fs1)<-rep("",(k[1])-1)
      fs2<-t(apply(bs1,1,sort))
      fs2<-delmin(fs2)
      rownames(fs2)<-paste("S",(i+1):(i+sk2), sep="")
      colnames(fs2)<-rep("",(k[2]-1))
      list(B1=list(fs1,fs2),B3=list(gs1,gs2),B4=A2)
    }
    #######################################################################
    # Obtaing set(s) of shifts by deleting smallest value of each group
    #######################################################################
    delmin<-function(z){
      fs<-c()
      n<-nrow(z)
      c<-ncol(z)-1
      for(i in 1:n){
        z1<-z[i,]
        z2<-z1[z1!=min(z1)]
        fs<-rbind(fs,z2)
      }
      return(fs)
    }
    ################################################################################
    # Selection of adjusted A and the set(s) of shifts to obtain Circular  
    # balance neighbor design for two different block size.
    ################################################################################
    # D=1: Circular Balanced Neighbor Designs
    # D=2: Circular Strongly Balanced Neighbor Designs
    #   K: Vector of two different block sizes
    #   i: Number of sets of shifts for K1
    # Sk2: Number of sets of shifts for K2
    CBND_2diffsize<-function(k,i,D=1,sk2=1){
      if(length(k)>2 | length(k)<2){stop("length(k)=2 ")}
      if(any(k<=2)!=0) stop("k=Block size: Each block size must be greater than 3")
      if(i<=0) stop("i= Must be a positive integer")
      if(k[1]<k[2]) stop("k1>K2")
      setClass( "stat_test", representation("list"))
      setMethod("show", "stat_test", function(object) {
        row <- paste(rep("=", 51), collapse = "")
        cat(row, "\n")
        
        row <- paste(rep("=", 51), collapse = "")
        cat(row, "\n")
        print(object$S[[1]])
        cat("\n")
        print(object$S[[2]])
      })
      if(D==2 & sk2==1){  
        v=2*i*k[1]+2*k[2]-1; m=(v-1)/2
        if(m%%8==0){
          j=m/8
          if(j<1) {return("Conditions are not satisfied for CSBND")}
          A=c(0:(j-1),(j+1):m,(v-j))
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
        if(m%%8==1){
          j=(m-1)/8
          if(j<1) {return("Conditions are not satisfied for CSBND")}
          A=c(0:(3*j),(3*j+2):(m-1),(m+1),(v-(3*j+1)))
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
        if(m%%8==2){
          j=(m-2)/8
          if(j<1) {return("Conditions are not satisfied for CSBND")}
          A=c(0:(5*j+1),(5*j+3):(m-1),(m+1),(v-(5*j+2)))
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
        if(m%%8==3){
          j=(m-3)/8
          if(j<0) {return("Conditions are not satisfied for CSBND")}
          A=c(0:(m-j-1),(m-j+1):m,(v-(m-j)))
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
        if(m%%8==4){
          j=(m-4)/8
          if(j<0) {return("Conditions are not satisfied for CSBND")}
          A=c(0:j,(j+2):(m-1),(m+1),(v-(j+1)))
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
        if(m%%8==5){
          j=(m-5)/8
          if(j<0) {return("Conditions are not satisfied for CSBND")}
          A=c(0:(3*j+1),(3*j+3):(m),(v-(3*j+2)))
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
        if(m%%8==6){
          j=(m-6)/8
          if(j<0) {return("Conditions are not satisfied for CSBND")}
          A=c(0:(5*j+3),(5*j+5):(m),(v-(5*j+4)))
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
        if(m%%8==7){
          j=(m-7)/8
          if(j<1) {return("Conditions are not satisfied for CSBND")}
          A=c(0:(m-j-1),(m-j+1):(m-1),(m+1),(v-(m-j)))
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        } 
      }
      if(D==2 & sk2==2){
        v=2*i*k[1]+4*k[2]-1; m=(v-1)/2
        if(m%%8==0){
          j=m/8
          if(j<1) {return("Conditions are not satisfied for CSBND")}
          A=c(0:(j-1),(j+1):m,(v-j))
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
        if(m%%8==1){
          j=(m-1)/8
          if(j<1) {return("Conditions are not satisfied for CSBND")}
          A=c(0:(3*j),(3*j+2):(m-1),(m+1),(v-(3*j+1)))
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
        if(m%%8==2){
          j=(m-2)/8
          if(j<1) {return("Conditions are not satisfied for CSBND")}
          A=c(0:(5*j+1),(5*j+3):(m-1),(m+1),(v-(5*j+2)))
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
        if(m%%8==3){
          j=(m-3)/8
          if(j<0) {return("Conditions are not satisfied for CSBND")}
          A=c(0:(m-j-1),(m-j+1):m,(v-(m-j)))
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
        if(m%%8==4){
          j=(m-4)/8
          if(j<0) {return("Conditions are not satisfied for CSBND")}
          A=c(0:j,(j+2):(m-1),(m+1),(v-(j+1)))
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
        if(m%%8==5){
          j=(m-5)/8
          if(j<0) {return("Conditions are not satisfied for CSBND")}
          A=c(0:(3*j+1),(3*j+3):(m),(v-(3*j+2)))
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
        if(m%%8==6){
          j=(m-6)/8
          if(j<0) {return("Conditions are not satisfied for CSBND")}
          A=c(0:(5*j+3),(5*j+5):(m),(v-(5*j+4)))
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
        if(m%%8==7){
          j=(m-7)/8
          if(j<1) {return("Conditions are not satisfied for CSBND")}
          A=c(0:(m-j-1),(m-j+1):(m-1),(m+1),(v-(m-j)))
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
      }
      if(D==1 & sk2==1){
        v=2*i*k[1]+2*k[2]+1; m=(v-1)/2
        if(m%%8==0){
          j=m/8
          if(j<1) {return("Conditions are not satisfied for CBNDs")}
          A=c(1:(j-1),(j+1):m,(v-j))
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
        if(m%%8==1){
          j=(m-1)/8
          if(j<1) {return("Conditions are not satisfied for CBNDs")}
          A=c(1:(3*j),(3*j+2):(m-1),(m+1),(v-(3*j+1)))
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
        if(m%%8==2){
          j=(m-2)/8
          if(j<1) {return("Conditions are not satisfied for CBNDs")}
          A=c(1:(5*j+1),(5*j+3):(m-1),(m+1),(v-(5*j+2)))
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
        if(m%%8==3){
          j=(m-3)/8
          if(j<0) {return("Conditions are not satisfied for CBNDs")}
          A=c(1:(m-j-1),(m-j+1):m,(v-(m-j)))
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
        if(m%%8==4){
          j=(m-4)/8
          if(j<0) {return("Conditions are not satisfied for CBNDs")}
          A=c(1:j,(j+2):(m-1),(m+1),(v-(j+1)))
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
        if(m%%8==5){
          j=(m-5)/8
          if(j<0) {return("Conditions are not satisfied for CBNDs")}
          A=c(1:(3*j+1),(3*j+3):(m),(v-(3*j+2)))
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
        if(m%%8==6){
          j=(m-6)/8
          if(j<0) {return("Conditions are not satisfied for CBNDs")}
          A=c(1:(5*j+3),(5*j+5):(m),(v-(5*j+4)))
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
        if(m%%8==7){
          j=(m-7)/8
          if(j<1) {return("Conditions are not satisfied for CBNDs")}
          A=c(1:(m-j-1),(m-j+1):(m-1),(m+1),(v-(m-j)))
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }   
      }   
      if(D==1 & sk2==2){
        v=2*i*k[1]+4*k[2]+1; m=(v-1)/2
        if(m%%8==0){
          j=m/8
          if(j<1) {return("Conditions are not satisfied for CBNDs")}
          A=c(1:(j-1),(j+1):m,(v-j))
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
        if(m%%8==1){
          j=(m-1)/8
          if(j<1) {return("Conditions are not satisfied for CBNDs")}
          A=c(1:(3*j),(3*j+2):(m-1),(m+1),(v-(3*j+1)))
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
        if(m%%8==2){
          j=(m-2)/8
          if(j<1) {return("Conditions are not satisfied for CBNDs")}
          A=c(1:(5*j+1),(5*j+3):(m-1),(m+1),(v-(5*j+2)))
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
        if(m%%8==3){
          j=(m-3)/8
          if(j<0) {return("Conditions are not satisfied for CBNDs")}
          A=c(1:(m-j-1),(m-j+1):m,(v-(m-j)))
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
        if(m%%8==4){
          j=(m-4)/8
          if(j<0) {return("Conditions are not satisfied for CBNDs")}
          A=c(1:j,(j+2):(m-1),(m+1),(v-(j+1)))
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
        if(m%%8==5){
          j=(m-5)/8
          if(j<0) {return("Conditions are not satisfied for CBNDs")}
          A=c(1:(3*j+1),(3*j+3):(m),(v-(3*j+2)))
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
        if(m%%8==6){
          j=(m-6)/8
          if(j<0) {return("Conditions are not satisfied for CBNDs")}
          A=c(1:(5*j+3),(5*j+5):(m),(v-(5*j+4)))
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
        if(m%%8==7){
          j=(m-7)/8
          if(j<1) {return("Conditions are not satisfied for CBNDs")}
          A=c(1:(m-j-1),(m-j+1):(m-1),(m+1),(v-(m-j)))
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }   
      }
      if(D==12 & sk2==1){  
        v=2*i*k[1]+2*k[2]-3; m=(v-1)/2
        if(m%%8==0){
          j=m/8
          if(j<1) {return("Conditions are not satisfied for CSBGND")}
          A=c(0:(j-1),(j+1),(j+2):(m-1),(m+1),(2*m-1-j),1)
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
        if(m%%8==1){
          j=(m-1)/8
          if(j<1) {return("Conditions are not satisfied for CSBGND")}
          A=c(0:(3*j),(3*j+2),(3*j+3):m,(2*m-3*j),1)
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
        if(m%%8==2){
          j=(m-2)/8
          if(j<1) {return("Conditions are not satisfied for CSBGND")}
          A=c(0:(5*j+1),(5*j+3),(5*j+4):m,(2*m-5*j-1),1)
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
        if(m%%8==3){
          j=(m-3)/8
          if(j<0) {return("Conditions are not satisfied for CSBGND")}
          A=c(0:(m-j),(m-j+2),(m-j+3):m,(m+j),1)
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
        if(m%%8==4){
          j=(m-4)/8
          if(j<0) {return("Conditions are not satisfied for CSBGND")}
          A=c(0:j,(j+2),(j+3):m,(2*m-j),1)
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
        if(m%%8==5){
          j=(m-5)/8
          if(j<0) {return("Conditions are not satisfied for CSBGND")}
          A=c(0:(3*j+2),(3*j+4),(3*j+5):(m-1),(m+1),(2*m-3*j-2),1)
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
        if(m%%8==6){
          j=(m-6)/8
          if(j<0) {return("Conditions are not satisfied for CSBGND")}
          A=c(0:(5*j+4),(5*j+6),(5*j+7):(m-1),(m+1),(2*m-5*j-4),1)
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
        if(m%%8==7){
          j=(m-7)/8
          if(j<1) {return("Conditions are not satisfied for CSBGND")}
          A=c(0:(m-j-1),(m-j+1),(m+j+2):m,(m+j+1),1)
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        } 
      }
      if(D==12 & sk2==2){
        v=2*i*k[1]+4*k[2]-3; m=(v-1)/2
        if(m%%8==0){
          j=m/8
          if(j<1) {return("Conditions are not satisfied for CSBGND")}
          A=c(0:(j-1),(j+1),(j+2):(m-1),(m+1),(2*m-1-j),1)
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
        if(m%%8==1){
          j=(m-1)/8
          if(j<1) {return("Conditions are not satisfied for CSBGND")}
          A=c(0:(3*j),(3*j+2),(3*j+3):m,(2*m-3*j),1)
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
        if(m%%8==2){
          j=(m-2)/8
          if(j<1) {return("Conditions are not satisfied for CSBGND")}
          A=c(0:(5*j+1),(5*j+3),(5*j+4):m,(2*m-5*j-1),1)
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
        if(m%%8==3){
          j=(m-3)/8
          if(j<0) {return("Conditions are not satisfied for CSBGND")}
          A=c(0:(m-j),(m-j+2),(m-j+3):m,(m+j),1)
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
        if(m%%8==4){
          j=(m-4)/8
          if(j<0) {return("Conditions are not satisfied for CSBGND")}
          A=c(0:j,(j+2),(j+3):m,(2*m-j),1)
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
        if(m%%8==5){
          j=(m-5)/8
          if(j<0) {return("Conditions are not satisfied for CSBGND")}
          A=c(0:(3*j+2),(3*j+4),(3*j+5):(m-1),(m+1),(2*m-3*j-2),1)
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
        if(m%%8==6){
          j=(m-6)/8
          if(j<0) {return("Conditions are not satisfied for CSBGND")}
          A=c(0:(5*j+4),(5*j+6),(5*j+7):(m-1),(m+1),(2*m-5*j-4),1)
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
        if(m%%8==7){
          j=(m-7)/8
          if(j<1) {return("Conditions are not satisfied for CSBGND")}
          A=c(0:(m-j-1),(m-j+1),(m+j+2):m,(m+j+1),1)
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
      }
      if(D==11 & sk2==1){
        v=2*i*k[1]+2*k[2]+3; m=(v-1)/2
        if(m%%8==0){
          j=m/8
          if(j<1) {return("Conditions are not satisfied for CPBND")}
          A=c(2:j,(j+2),(j+3):(m-2),m,(m+2),(2*m-j))
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
        if(m%%8==1){
          j=(m-1)/8
          if(j<1) {return("Conditions are not satisfied for CPBND")}
          A=c(2:(3*j-1),(3*j+1),(3*j+2):m,(2*m-(3*j-1)))
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
        if(m%%8==2){
          j=(m-2)/8
          if(j<1) {return("Conditions are not satisfied for CPBND")}
          A=c(2:(5*j+2),(5*j+3):m,(2*m-(5*j)))
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
        if(m%%8==3){
          j=(m-3)/8
          if(j<0) {return("Conditions are not satisfied for CPBND")}
          A=c(2:(m-j-1),(m-j+1),(m-j+2):(m+1),(m+2),(m+j+1))
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
        if(m%%8==4){
          j=(m-4)/8
          if(j<0) {return("Conditions are not satisfied for CPBND")}
          A=c(2:(j+1),(j+3),(j+4):(m-2),(m+1),(m+2),(2*m-j-1))
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
        if(m%%8==5){
          j=(m-5)/8
          if(j<0) {return("Conditions are not satisfied for CPBND")}
          A=c(2:(3*j+1),(3*j+3),(3*j+4):(m-1),(m+1),(2*m-3*j-1))
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
        if(m%%8==6){
          j=(m-6)/8
          if(j<0) {return("Conditions are not satisfied for CPBND")}
          A=c(2:(5*j+3),(5*j+5),(5*j+6):(m-1),(m+1),(2*m-5*j-3))
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
        if(m%%8==7){
          j=(m-7)/8
          if(j<1) {return("Conditions are not satisfied for CPBND")}
          A=c(2:(m-j-2),(m-j),(m-j+1):m,(m+j+2))
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }   
      }   
      if(D==11 & sk2==2){
        v=2*i*k[1]+4*k[2]+3; m=(v-1)/2
        if(m%%8==0){
          j=m/8
          if(j<1) {return("Conditions are not satisfied for CPBND")}
          A=c(2:j,(j+2),(j+3):(m-2),m,(m+2),(2*m-j))
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
        if(m%%8==1){
          j=(m-1)/8
          if(j<1) {return("Conditions are not satisfied for CPBND")}
          A=c(2:(3*j-1),(3*j+1),(3*j+2):m,(2*m-(3*j-1)))
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
        if(m%%8==2){
          j=(m-2)/8
          if(j<1) {return("Conditions are not satisfied for CPBND")}
          A=c(2:(5*j+2),(5*j+3):m,(2*m-(5*j)))
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
        if(m%%8==3){
          j=(m-3)/8
          if(j<0) {return("Conditions are not satisfied for CPBND")}
          A=c(2:(m-j-1),(m-j+1),(m-j+2):(m+1),(m+2),(m+j+1))
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
        if(m%%8==4){
          j=(m-4)/8
          if(j<0) {return("Conditions are not satisfied for CPBND")}
          A=c(2:(j+1),(j+3),(j+4):(m-2),(m+1),(m+2),(2*m-j-1))
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
        if(m%%8==5){
          j=(m-5)/8
          if(j<0) {return("Conditions are not satisfied for CPBND")}
          A=c(2:(3*j+1),(3*j+3),(3*j+4):(m-1),(m+1),(2*m-3*j-1))
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
        if(m%%8==6){
          j=(m-6)/8
          if(j<0) {return("Conditions are not satisfied for CPBND")}
          A=c(2:(5*j+3),(5*j+5),(5*j+6):(m-1),(m+1),(2*m-5*j-3))
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
        if(m%%8==7){
          j=(m-7)/8
          if(j<1) {return("Conditions are not satisfied for CPBND")}
          A=c(2:(m-j-2),(m-j),(m-j+1):m,(m+j+2))
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }   
      }
      
      if(D==5 & sk2==1){  
        v=2*i*k[1]+2*k[2]; m=(v-2)/2
        if(m%%4==2){
          A=1:(m+1)
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
        if(m%%4==3){
          A<-c(1:((m-3)/4),((m+5)/4),((m+9)/4):m,(m+1),(7*(m+1)/4))
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
        if(m%%4==0 | m%%4==1){return("The minimal CWBNDs in which v/2 unordered pair apears twice cannot be constructed for v=2ik1+2k2 and k1= block size and k2= block size")}
      }
      if(D==5 & sk2==2){
        v=2*i*k[1]+4*k[2]; m=(v-2)/2
        if(m%%4==2){
          A=1:(m+1)
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
        if(m%%4==3){
          A<-c(1:((m-3)/4),((m+5)/4),((m+9)/4):m,(m+1),(7*(m+1)/4))
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
        if(m%%4==0 | m%%4==1){return("The minimal CWBNDs in which v/2 unordered pair apears twice cannot be constructed for v=2ik1+2k2 and k1=block size and k2= block size")}
      }
      if(D==6 & sk2==1){
        v=2*i*k[1]+2*k[2]-2; m=(v-2)/2
        if(m%%4==0){
          A=c(1:(m-1),(m-1),(m+1),(m+2))
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
        if(m%%4==1){
          A=c(1:((m-5)/4),((m+3)/4),((m+7)/4):m,(m+1),m,((7*m+9)/4))
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
        if(m%%4==2){
          A=c(1:((m-2)/2),((m+2)/2),((m+4)/2):m,(m+1),m,((3*m+4)/2))
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
        if(m%%4==3){
          A=c(1:((3*m-5)/4),((3*m+3)/4),((3*m+7)/4):(m-1), m,(m+1),(m-1), ((5*m+9)/4))
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
      }  
      if(D==6 & sk2==2){
        v=2*i*k[1]+4*k[2]-2; m=(v-2)/2
        if(m%%4==0){
          A=c(1:(m-1),(m-1),(m+1),(m+2))
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
        if(m%%4==1){
          A=c(1:((m-5)/4),((m+3)/4),((m+7)/4):m,(m+1),m,((7*m+9)/4))
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
        if(m%%4==2){
          A=c(1:((m-2)/2),((m+2)/2),((m+4)/2):m,(m+1),m,((3*m+4)/2))
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
        if(m%%4==3){
          A=c(1:((3*m-5)/4),((3*m+3)/4),((3*m+7)/4):(m-1), m,(m+1),(m-1), ((5*m+9)/4))
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
      }
      if(D==7 & sk2==1){
        v=2*i*k[1]+2*k[2] ; m=(v-2)/2
        if(m%%4==0){
          A<-0:m
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
        if(m%%4==3){
          A<-c(0:((3*m-1)/4),((3*m+7)/4),((3*m+11)/4):m, (5*(m+1)/4))
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
        if(m%%4==1 | m%%4==2){return("The minimal CSPBNDs in which v/2 unordered pair cannot be constructed for V=2ik1+2k2+2 and k1= block size and k2= block size")}
      } 
      if(D==7 & sk2==2){
        v=2*i*k[1]+4*k[2] ; m=(v-2)/2
        if(m%%4==0){
          A<-0:m
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
        if(m%%4==3){
          A<-c(0:((3*m-1)/4),((3*m+7)/4),((3*m+11)/4):m, (5*(m+1)/4))
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
        if(m%%4==1 | m%%4==2){return("The minimal CSPBNDs in which v/2 unordered pair cannot be constructed for V=2ik1+2k2+2 and k1= block size and k2= block size")}
      }
      if(D==8 & sk2==1){
        v= 2*i*k[1]+2*k[2]+2 ; m=(v-2)/2 
        if(m%%4==0){
          A=c(0:(m/2),((m+4)/2),((m+6)/2):(m-1),((3*m+2)/2))
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
        if(m%%4==1){
          A=c(0:((3*m+1)/4),((3*m+9)/4),((3*m+13)/4):(m-1),((5*m+3)/4))
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
        if(m%%4==2){
          A=c(0,2:(m-2),m,(2*m+1))
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
        if(m%%4==3){
          A=c(0:((m+1)/4),((m+9)/4),((m+13)/4):(m-2), m, ((7*m+3)/4))
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
      } 
      if(D==8 & sk2==2){
        v= 2*i*k[1]+4*k[2]+2 ; m=(v-2)/2 
        if(m%%4==0){
          A=c(0:(m/2),((m+4)/2),((m+6)/2):(m-1),((3*m+2)/2))
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
        if(m%%4==1){
          A=c(0:((3*m+1)/4),((3*m+9)/4),((3*m+13)/4):(m-1),((5*m+3)/4))
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
        if(m%%4==2){
          A=c(0,2:(m-2),m,(2*m+1))
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
        if(m%%4==3){
          A=c(0:((m+1)/4),((m+9)/4),((m+13)/4):(m-2), m, ((7*m+3)/4))
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
      }
      if(D==9 & sk2==1){  
        v=2*i*k[1]+2*k[2]-2; m=(v-2)/2
        if(m%%4==2){
          A=0:(m+1)
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
        if(m%%4==3){
          A<-c(0:((m-3)/4),((m+5)/4),((m+9)/4):m,(m+1),(7*(m+1)/4))
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
        if(m%%4==0 | m%%4==1){return("The minimal CSBGNDs in which v/2 unordered pair apears twice cannot be constructed for v=2ik1+2k2 and k1= block size and k2= block size")}
      }  
      if(D==9 & sk2==2){
        v=2*i*k[1]+4*k[2]-2; m=(v-2)/2
        if(m%%4==2){
          A=0:(m+1)
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
        if(m%%4==3){
          A<-c(0:((m-3)/4),((m+5)/4),((m+9)/4):m,(m+1),(7*(m+1)/4))
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
        if(m%%4==0 | m%%4==1){return("The minimal CSBGNDs in which v/2 unordered pair apears twice cannot be constructed for v=2ik1+2k2 and k1=block size and k2= block size")}
      }
      if(D==10 & sk2==1){
        v=2*i*k[1]+2*k[2]-4; m=(v-2)/2
        if(m%%4==0){
          A=c(0:(m-1),(m-1),(m+1),(m+2))
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
        if(m%%4==1){
          A=c(0:((m-5)/4),((m+3)/4),((m+7)/4):m,(m+1),m,((7*m+9)/4))
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
        if(m%%4==2){
          A=c(0:((m-2)/2),((m+2)/2),((m+4)/2):m,(m+1),m,((3*m+4)/2))
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
        if(m%%4==3){
          A=c(0:((3*m-5)/4),((3*m+3)/4),((3*m+7)/4):(m-1), m,(m+1),(m-1), ((5*m+9)/4))
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
      }  
      if(D==10 & sk2==2){
        v=2*i*k[1]+4*k[2]-4; m=(v-2)/2
        if(m%%4==0){
          A=c(0:(m-1),(m-1),(m+1),(m+2))
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
        if(m%%4==1){
          A=c(0:((m-5)/4),((m+3)/4),((m+7)/4):m,(m+1),m,((7*m+9)/4))
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
        if(m%%4==2){
          A=c(0:((m-2)/2),((m+2)/2),((m+4)/2):m,(m+1),m,((3*m+4)/2))
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
        if(m%%4==3){
          A=c(0:((3*m-5)/4),((3*m+3)/4),((3*m+7)/4):(m-1), m,(m+1),(m-1), ((5*m+9)/4))
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
      }
      if(D==3 & sk2==1){  
        v=2*i*k[1]+2*k[2]+2 ; m=(v-2)/2
        if(m%%4==0){
          A<-c(1:m)
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
        if(m%%4==3){
          A<-c(1:((3*m-1)/4),((3*m+7)/4),((3*m+11)/4):m, (5*(m+1)/4))
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
        if(m%%4==1 | m%%4==2){return("The minimal CPNDs in which v/2 unordered pair cannot be constructed for V=2ik1+2k2+2 and k1= block size and k2= block size")}
      }  
      if(D==3 & sk2==2){
        v=2*i*k[1]+4*k[2]+2 ; m=(v-2)/2
        if(m%%4==0){
          A<-c(1:m)
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
        if(m%%4==3){
          A<-c(1:((3*m-1)/4),((3*m+7)/4),((3*m+11)/4):m, (5*(m+1)/4))
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
        if(m%%4==1 | m%%4==2){return("The minimal CPNDs in which v/2 unordered pair cannot be constructed for V=2ik1+2k2+2 and k1= block size and k2= block size")}
      }  
      if(D==4 & sk2==1){
        v= 2*i*k[1]+2*k[2]+4 ; m=(v-2)/2 
        if(m%%4==0){
          A=c(1:(m/2),((m+4)/2),((m+6)/2):(m-1),((3*m+2)/2))
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
        if(m%%4==1){
          A=c(1:((3*m+1)/4),((3*m+9)/4),((3*m+13)/4):(m-1),((5*m+3)/4))
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
        if(m%%4==2){
          A=c(2:(m-2),m,(2*m+1))
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
        if(m%%4==3){
          A=c(1:((m+1)/4),((m+9)/4),((m+13)/4):(m-2), m, ((7*m+3)/4))
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
      }
      if(D==4 & sk2==2){
        v= 2*i*k[1]+4*k[2]+4  ; m=(v-2)/2 
        if(m%%4==0){
          A=c(1:(m/2),((m+4)/2),((m+6)/2):(m-1),((3*m+2)/2))
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
        if(m%%4==1){
          A=c(1:((3*m+1)/4),((3*m+9)/4),((3*m+13)/4):(m-1),((5*m+3)/4))
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
        if(m%%4==2){
          A=c(2:(m-2),m,(2*m+1))
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
        if(m%%4==3){
          A=c(1:((m+1)/4),((m+9)/4),((m+13)/4):(m-2), m, ((7*m+3)/4))
          A1<-grouping2(A,k,v,i,sk2)
          A2<-c(v,k);names(A2)<-c("V","K1","K2")
          x<-list(S=A1$B1,G=A1$B3,R=A2,A=A)
        }
      }
      new("stat_test", x) 
    }
    
    design_CBND<-function(H){
      setClass( "CBND_design", representation("list"))
      setMethod("show", "CBND_design", function(object) {
        row <- paste(rep("=", 51), collapse = "")
        cat(row, "\n")
        cat("Following is minimal CBND for", "v=" ,object$R[1], "and","k=",object$R[2], "\n")
        row <- paste(rep("=", 51), collapse = "")
        cat(row, "\n")
        for(i in 1:length(ss)){
          W<-ss[[i]]
          nr<-dim(W)[1]
          for(j in 1:nr){
            print(object$Design[[i]][[j]])
            cat("\n\n")
          }}
      }) 
      v<-H$R[1]
      k<-H$R[2]
      ss<-H$S  
      treat<-(1:v)-1
      fn<-(1:v)
      G<-list()
      for(j in 1:length(ss)){ 
        W<-ss[[j]]
        nr<-dim(W)[1]
        nc<-dim(W)[2]
        D<-list()
        for(i in 1:nr){
          dd<-c()
          d1<-matrix(treat,(nc+1),v,byrow = T)
          ss1<-cumsum(c(0,W[i,]))
          dd2<-d1+ss1
          dd<-rbind(dd,dd2)
          rr<-dd[which(dd>=v)]%%v
          dd[which(dd>=v)]<-rr
          colnames(dd)<-paste("B",fn, sep="")
          rownames(dd)<-rep("",(nc+1))
          fn<-fn+v
          D[[i]]<-dd
        }
        G[[j]]<-D
      }
      x<-list(Design=G,R=H$R)
      new("CBND_design", x)
    }
    Create_Designs<-function(v,k1,k2,i){
      if(v%%2!=0){
        if((v-1-(2*k2))%%(2*k1)==0 &(i=(v-1-(2*k2))/(2*k1)) &(i%%1==0)&(i>0)){
          cat("MCBNDs is possible with one set of k2","v=",v,"k1=",k1,"k2=",k2,"i=",i)
          print(H<-CBND_2diffsize(k=c(k1,k2),i=i,D=1,sk2=1))
        }
          if((v-1-(4*k2))%%(2*k1)==0 &(i=(v-1-(4*k2))/(2*k1)) & (i%%1==0) &(i>0) ){
            cat("MCBNDs is possible with two sets of k2","v=",v,"k1=",k1,"k2=",k2,"i=",i)
            print(H<-CBND_2diffsize(k=c(k1,k2),i=i,D=1,sk2=2))
          }
        if(((v+1-(2*k2))%%(2*k1)==0) &(i=(v+1-(2*k2))/(2*k1)) & (i%%1==0) &(i>0)){
                cat("MCSBNDs is possible with one set of shift k2","v=",v,"k1=",k1,"k2=",k2,"i=",i)
                print(H<-CBND_2diffsize(k=c(k1,k2),i=i,D=2,sk2=1))
              }else
                if((v+1-(4*k2))%%(2*k1)==0&(i=(v+1-(4*k2))/(2*k1)) & (i%%1==0) &(i>0)){
                  cat("MCSBND is possible with two     set of k2","v=",v,"k1=",k1,"k2=",k2,"i=",i)
                  print(H<-CBND_2diffsize(k=c(k1,k2),i=i,D=2,sk2=2))
                }
                  if((v-3-(2*k2))%%(2*k1)==0 &(i=(v-3-(2*k2))/(2*k1)) &(i%%1==0)&(i>0)){
                    cat("MCPBND-Odd is possible with one set of k2","v=",v,"k1=",k1,"k2=",k2,"i=",i)
                    print(H<-CBND_2diffsize(k=c(k1,k2),i=i,D=11,sk2=1))
                  }else
                    if((v-3-(4*k2))%%(2*k1)==0 &(i=(v-3-(4*k2))/(2*k1)) & (i%%1==0) &(i>0) ){
                      cat("MCPBND-Odd is possible with two sets of k2","v=",v,"k1=",k1,"k2=",k2,"i=",i)
                      print(H<-CBND_2diffsize(k=c(k1,k2),i=i,D=11,sk2=2))
                    }
                      if(((v+3-(2*k2))%%(2*k1)==0) &(i=(v+3-(2*k2))/(2*k1)) & (i%%1==0) &(i>0)){
                        cat("MCSBGND-Odd is possible with one set of shift k2","v=",v,"k1=",k1,"k2=",k2,"i=",i)
                        print(H<-CBND_2diffsize(k=c(k1,k2),i=i,D=12,sk2=1))
                      }else
                        if((v+3-(4*k2))%%(2*k1)==0&(i=(v+3-(4*k2))/(2*k1)) & (i%%1==0) &(i>0)){
                          cat("MCSBGND-Odd is possible with two     set of k2","v=",v,"k1=",k1,"k2=",k2,"i=",i)
                          print(H<-CBND_2diffsize(k=c(k1,k2),i=i,D=12,sk2=2))
                        }
        if((i%%1<0) |(i%%1>0)){cat("Other designs are not possible","v=",v,"k1=",k1,"k2=",k2,"i=",i)
        }
          }else
      { if((v-(2*k2)-2)%%(2*k1)==0 &(i=(v-(2*k2)-2)/(2*k1)) & (i%%1==0) &(i>0) & (((v-2)/2)%%4==0 |((v-2)/2)%%4==3)){
        cat("MCPBNDs-I is possible with one set of k2","v=",v,"k1=",k1,"k2=",k2,"i=",i,rep("",15))
        print(H<-CBND_2diffsize(k=c(k1,k2),i=i,D=3,sk2=1))
      }else
        if((v-(4*k2)-2)%%(2*k1)==0&(i=(v-(4*k2)-2)/(2*k1)) & (i%%1==0) &(i>0) & (((v-2)/2)%%4==0 |((v-2)/2)%%4==3)){
          cat("MCPBNDs-I is possible with two set of k2","v=",v,"k1=",k1,"k2=",k2,"i=",i,rep("",15))
          print(H<-CBND_2diffsize(k=c(k1,k2),i=i,D=3,sk2=2))
        }
          if((v-(2*k2)-2)%%(2*k1)==0&(i=(v-(2*k2)-2)/(2*k1)) & (i%%1==0) &(i>0)){
            cat("MCSPBNDS-II is possible with one set of k2","v=",v,"k1=",k1,"k2=",k2,"i=",i,rep("",15))
            print(H<-CBND_2diffsize(k=c(k1,k2),i=i,D=8,sk2=1))
          }else
            if((v-(4*k2)-2)%%(2*k1)==0 & (i=(v-(4*k2)-2)/(2*k1)) & (i%%1==0) &(i>0)){
              cat(" MCSPBNDS-II is possible with two set of k2","v=",v,"k1=",k1,"k2=",k2,"i=",i,rep("",15))
              print(H<-CBND_2diffsize(k=c(k1,k2),i=i,D=8,sk2=2))
            }
              if((v-(2*k2)-4)%%(2*k1)==0&(i=(v-(2*k2)-4)/(2*k1)) & (i%%1==0) &(i>0)){
                cat(" MCPBNDS-II is possible with one set of k2","v=",v,"k1=",k1,"k2=",k2,"i=",i,rep("",15))
                print(H<-CBND_2diffsize(k=c(k1,k2),i=i,D=4,sk2=1))
              }else
                if((v-(4*k2)-4)%%(2*k1)==0 & (i=(v-(4*k2)-4)/(2*k1)) & (i%%1==0) &(i>0)){
                  cat(" MCPBNDS-II is possible with two set of k2","v=",v,"k1=",k1,"k2=",k2,"i=",i,rep("",15))
                  print(H<-CBND_2diffsize(k=c(k1,k2),i=i,D=4,sk2=2))
                }
                  if((v-(2*k2))%%(2*k1)==0&(i=(v-(2*k2))/(2*k1)) & (i%%1==0) &(i>0)& (((v-2)/2)%%4==2 |((v-2)/2)%%4==3)) {
                    cat("MCWBNDs-I is possible with one set of k2", "v=" ,v,"k1=", k1,"k2=", k2, "i=",i,rep("  ",20))
                    print(H<-CBND_2diffsize(k=c(k1,k2),i=i,D=5,sk2=1))
                  }else
                    if((v-(4*k2))%%(2*k1)==0&(i=(v-(4*k2))/(2*k1)) & (i%%1==0) &(i>0)& (((v-2)/2)%%4==2 |((v-2)/2)%%4==3)) {
                      cat("MCWBNDs-I is possible with two set of k2", "v=" ,v,"k1=", k1,"k2=", k2, "i=",i,rep("  ",20))
                      print(H<-CBND_2diffsize(k=c(k1,k2),i=i,D=5,sk2=2))
                    }
                      if((v-(2*k2))%%(2*k1)==0&(i=(v-(2*k2))/(2*k1)) & (i%%1==0) &(i>0)& (((v-2)/2)%%4==0 |((v-2)/2)%%4==3)){
                        cat("MCSPBNDs-I is possible with one set of k2","v=",v,"k1=",k1,"k2=",k2,"i=",i,rep("",15))
                        print(H<-CBND_2diffsize(k=c(k1,k2),i=i,D=7,sk2=1))
                      }else
                        if((v-(4*k2))%%(2*k1)==0&(i=(v-(4*k2))/(2*k1)) & (i%%1==0) &(i>0)& (((v-2)/2)%%4==0 |((v-2)/2)%%4==3)){
                          cat("MCSPBNDs-I is possible with two set of k2","v=",v,"k1=",k1,"k2=",k2,"i=",i,rep("",15))
                          print(H<-CBND_2diffsize(k=c(k1,k2),i=i,D=7,sk2=2))
                        }
                              if((v-(2*k2)+2)%%(2*k1)==0&(i=(v-(2*k2)+2)/(2*k1)) & (((v-2)/2)%%4==2 |((v-2)/2)%%4==3)& (i%%1==0) &(i>0)){
                                cat("MCSGBNDS-I is possible with one set of k2","v=",v,"k1=",k1,"k2=",k2,"i=",i,rep("",15))
                                print(H<-CBND_2diffsize(k=c(k1,k2),i=i,D=9,sk2=1))
                              }else
                                if((v-(4*k2)+2)%%(2*k1)==0&(i=(v-(4*k2)+2)/(2*k1)) & (((v-2)/2)%%4==2 |((v-2)/2)%%4==3)& (i%%1==0) &(i>0)){
                                  cat("MCSGBNDS-I is possible with two set of k2","v=",v,"k1=",k1,"k2=",k2,"i=",i,rep("",15))
                                  print(H<-CBND_2diffsize(k=c(k1,k2),i=i,D=9,sk2=2))
                                } 
                                  if((v-(2*k2)+2)%%(2*k1)==0&(i=(v-(2*k2)+2)/(2*k1)) & (i%%1==0) &(i>0) ){
                                    cat(" MCWBNDs-II is possible with one set of k2","v=",v,"k1=",k1,"k2=",k2,"i=",i,rep("",15))
                                    print(H<-CBND_2diffsize(k=c(k1,k2),i=i,D=6,sk2=1))
                                  }else 
                                    if((v-(4*k2)+2)%%(2*k1)==0 & (i=(v-(4*k2)+2)/(2*k1)) & (i%%1==0) & (i>0)){
                                      cat("MCWBNDs-II is possible with two set of k2","v=",v,"k1=",k1,"k2=",k2,"i=",i,rep("",15))
                                      print(H<-CBND_2diffsize(k=c(k1,k2),i=i,D=6,sk2=2))
                                    }
                                      if((v-(2*k2)+4)%%(2*k1)==0&(i=(v-(2*k2)+4)/(2*k1)) & (i%%1==0) &(i>0)){
                                        cat(" MCSGBNDS-II is possible with one set of k2","v=",v,"k1=",k1,"k2=",k2,"i=",i,rep("",15))
                                        print(H<-CBND_2diffsize(k=c(k1,k2),i=i,D=10,sk2=1))
                                      }else
                                        if((v-(4*k2)+4)%%(2*k1)==0 & (i=(v-(4*k2)+4)/(2*k1)) & (i%%1==0) &(i>0)){
                                          cat(" MCSGBNDS-II is possible with two set of k2","v=",v,"k1=",k1,"k2=",k2,"i=",i,rep("",15))
                                          print(H<-CBND_2diffsize(k=c(k1,k2),i=i,D=10,sk2=2))
                                        }
                                      
                                          if((i%%1<0) |(i%%1>0)) cat("Other designs are not possible","v=",v,"k1=",k1,"k2=",k2,"i=",i," ")
      }
    }
    #MCBNDs and MCSBNDs, say here D1 & D2 respectively, for v odd
    
    (H<-Create_Designs(v=v, k1=k1, k2=k2))
    
  }else
    if(v>0 &k1>0&k2>0 &k3>0){################################################################################
      # CBND_3diffsize: Circular balance neighbor design for block of three different 
      # sizes (K1,k2 and k3)
      # Algorithm from paper:
      # Akbar Firdos,Mahmood Ul Hassan,Farrukh Jamal,Hurria Ali,Khadija Noreen and Rashid Ahmed Algorithms to Construct Minimal Circular Strongly. 
      # Coded by Fardos Akbar et al., 2023
      ################################################################################
      # Selection of i group of size K1 from adjusted A. The set of remaining 
      # (Unselected) elements are saved in object B2. 
      ################################################################################
      grouping1<-function(A,k,v,i){
        bs<-c()
        z=0;f=1
        A1=A
        while(f<=i){
          for(y in 1:5000){
            comp<-sample(1:length(A1),k)
            com<-A1[comp]
            cs<-sum(com)
            if(cs%%v==0){
              bs<-rbind(bs,com)
              A1<-A1[-comp]
              z<-z+1
              f=f+1
            }
            if(z==i) break
          }
          if(z<i) {bs<-c();z=0;f=1;A1=A}  
        }
        list(B1=bs,B2=A1)
      }
      ################################################################################
      # Selection of i group of size K1 from adjusted A and selection of required 
      # number of groups of size K2 from B2. The set of remaining (Unselected) 
      # elements are saved in B3.
      ################################################################################
      grouping2<-function(A,k,v,i,sk2){
        bs1<-c()
        j=i+sk2
        z=0;f=1
        A1=A
        while(f<=j){
          s<-grouping1(A1,k[1],v,i)
          A2<-s$B2
          z=i;f=f+i
          for(y in 1:2000){
            comp<-sample(1:length(A2),k[2])
            com<-A2[comp]
            cs<-sum(com)
            if(cs%%v==0){
              bs1<-rbind(bs1,com)
              A2<-A2[-comp]
              z<-z+1
              f=f+1
            }
            if(z==j) break
          }
          if(z<j) {bs1<-c();z=0;f=1;A1=A}  
        }
        list(B1=s$B1,B2=bs1,B3=A2)
      }
      ################################################################################
      # Selection of i group of size K1 from adjusted A, selection of required number 
      # of groups of size K2 from B2 and division of required number of groups of size
      # K3 from B3.
      ################################################################################
      grouping3<-function(A,k,v,i,sk2,sk3){
        bs1<-c()
        j=i+sk2+sk3
        z=0;f=1
        A1=A
        while(f<=j){
          s<-grouping2(A1,k,v,i,sk2)
          A3<-s$B3
          z=i+sk2;f=f+i+sk2
          for(y in 1:1000){
            comp<-sample(1:length(A3),k[3])
            com<-A3[comp]
            cs<-sum(com)
            if(cs%%v==0){
              bs1<-rbind(bs1,com)
              A3<-A3[-com]
              z<-z+1
              f=f+1
            }
            if(z==j) break
          }
          if(z<j) {bs1<-c();z=0;f=1;A1=A}  
        }
        gs1<-t(apply(s$B1,1,sort))
        gs1<-cbind(gs1,rowSums(gs1),rowSums(gs1)/v)
        rownames(gs1)<-paste("G",1:i, sep="")
        colnames(gs1)<-c(paste(1:k[1], sep=""),"sum" ,"sum/v")
        gs2<-t(apply(s$B2,1,sort))
        gs2<-cbind(gs2,rowSums(gs2),rowSums(gs2)/v)
        rownames(gs2)<-paste("G",(i+1):(i+sk2), sep="")
        colnames(gs2)<-c(paste(1:k[2], sep=""),"sum" ,"sum/v")
        gs3<-t(apply(bs1,1,sort))
        gs3<-cbind(gs3,rowSums(gs3),rowSums(gs3)/v)
        rownames(gs3)<-paste("G",(i+sk2+1):(i+sk2+sk3), sep="")
        colnames(gs3)<-c(paste(1:k[3], sep=""),"sum" ,"sum/v")
        fs1<-t(apply(s$B1,1,sort))
        fs1<-delmin(fs1)
        rownames(fs1)<-paste("S",1:i, sep="")
        colnames(fs1)<-rep("",(k[1])-1)
        fs2<-t(apply(s$B2,1,sort))
        fs2<-delmin(fs2)
        rownames(fs2)<-paste("S",(i+1):(i+sk2), sep="")
        colnames(fs2)<-rep("",(k[2])-1)
        fs3<-t(apply(bs1,1,sort))
        fs3<-delmin(fs3)
        rownames(fs3)<-paste("S",(i+sk2+1):(i+sk2+sk3), sep="")
        colnames(fs3)<-rep("",(k[3]-1))
        list(B1=list(fs1,fs2,fs3),B4=list(gs1,gs2,gs3),B5=A3)
      }
      #######################################################################
      # Obtain set(s) of shifts by deleting smallest value of each group
      #######################################################################
      delmin<-function(z){
        fs<-c()
        n<-nrow(z)
        c<-ncol(z)-1
        for(i in 1:n){
          z1<-z[i,]
          z2<-z1[z1!=min(z1)]
          fs<-rbind(fs,z2)
        }
        return(fs)
      }
      ################################################################################
      # Selection of adjusted A and the set(s) of shifts to obtain Circular  
      # balance neighbor design for three different block size.
      ################################################################################
      CBND_3diffsize<-function(k,i,D,sk2,sk3){
        if(length(k)>3 | length(k)<3){stop("length(k)=3")}
        if(any(k<=2)!=0) stop("k=Block size: Each block size must be greater than 3")
        if(i<=0) stop("i=must be a positive integer")
        if(k[1]<k[2] | k[2]<k[3] |  k[1]<k[3]  ) stop("k1>K2>K3")
        
        setClass( "stat_test", representation("list"))
        setMethod("show", "stat_test", function(object) {
          row <- paste(rep("", 52), collapse = "")
          cat(row, "\n")
          row <- paste(rep("=", 52), collapse = "")
          cat(row, "\n")
          
          print(object$S[[1]])
          print(object$S[[2]])
          print(object$S[[3]])
        })
        if(D==2 & sk2==1 & sk3==1){  
          v=2*i*k[1]+2*k[2]+2*k[3]-1 ; m=(v-1)/2
          if(m%%8==0){
            j=m/8
            if(j<1) {return("Conditions are not satisfied for CSBND")}
            A=c(0:(j-1),(j+1):m,(v-j))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==1){
            j=(m-1)/8
            if(j<1) {return("Conditions are not satisfied for CSBND")}
            A=c(0:(3*j),(3*j+2):(m-1),(m+1),(v-(3*j+1)))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==2){
            j=(m-2)/8
            if(j<1) {return("Conditions are not satisfied for CSBND")}
            A=c(0:(5*j+1),(5*j+3):(m-1),(m+1),(v-(5*j+2)))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==3){
            j=(m-3)/8
            if(j<0) {return("Conditions are not satisfied for CSBND")}
            A=c(0:(m-j-1),(m-j+1):m,(v-(m-j)))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==4){
            j=(m-4)/8
            if(j<0) {return("Conditions are not satisfied for CSBND")}
            A=c(0:j,(j+2):(m-1),(m+1),(v-(j+1)))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==5){
            j=(m-5)/8
            if(j<0) {return("Conditions are not satisfied for CSBND")}
            A=c(0:(3*j+1),(3*j+3):(m),(v-(3*j+2)))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==6){
            j=(m-6)/8
            if(j<0) {return("Conditions are not satisfied for CSBND")}
            A=c(0:(5*j+3),(5*j+5):(m),(v-(5*j+4)))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==7){
            j=(m-7)/8
            if(j<1) {return("Conditions are not satisfied for CSBND")}
            A=c(0:(m-j-1),(m-j+1):(m-1),(m+1),(v-(m-j)))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
        }
        if(D==1 & sk2==1 & sk3==1){
          v= 2*i*k[1]+2*k[2]+2*k[3]+1; m=(v-1)/2
          if(m%%8==0){
            j=m/8
            if(j<1) {return("Conditions are not satisfied for CBNDs")}
            A=c(1:(j-1),(j+1):m,(v-j))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==1){
            j=(m-1)/8
            if(j<1) {return("Conditions are not satisfied for CBNDs")}
            A=c(1:(3*j),(3*j+2):(m-1),(m+1),(v-(3*j+1)))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==2){
            j=(m-2)/8
            if(j<1) {return("Conditions are not satisfied for CBNDs")}
            A=c(1:(5*j+1),(5*j+3):(m-1),(m+1),(v-(5*j+2)))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==3){
            j=(m-3)/8
            if(j<0) {return("Conditions are not satisfied for CBNDs")}
            A=c(1:(m-j-1),(m-j+1):m,(v-(m-j)))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==4){
            j=(m-4)/8
            if(j<0) {return("Conditions are not satisfied for CBNDs")}
            A=c(1:j,(j+2):(m-1),(m+1),(v-(j+1)))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==5){
            j=(m-5)/8
            if(j<0) {return("Conditions are not satisfied for CBNDs")}
            A=c(1:(3*j+1),(3*j+3):(m),(v-(3*j+2)))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==6){
            j=(m-6)/8
            if(j<0) {return("Conditions are not satisfied for CBNDs")}
            A=c(1:(5*j+3),(5*j+5):(m),(v-(5*j+4)))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==7){
            j=(m-7)/8
            if(j<1) {return("Conditions are not satisfied for CBNDs")}
            A=c(1:(m-j-1),(m-j+1):(m-1),(m+1),(v-(m-j)))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
        }  
        if(D==2 & sk2==2 & sk3==2){  
          v=2*i*k[1]+4*k[2]+4*k[3]-1 ; m=(v-1)/2
          if(m%%8==0){
            j=m/8
            if(j<1) {return("Conditions are not satisfied for CSBND")}
            A=c(0:(j-1),(j+1):m,(v-j))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==1){
            j=(m-1)/8
            if(j<1) {return("Conditions are not satisfied for CSBND")}
            A=c(0:(3*j),(3*j+2):(m-1),(m+1),(v-(3*j+1)))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==2){
            j=(m-2)/8
            if(j<1) {return("Conditions are not satisfied for CSBND")}
            A=c(0:(5*j+1),(5*j+3):(m-1),(m+1),(v-(5*j+2)))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==3){
            j=(m-3)/8
            if(j<0) {return("Conditions are not satisfied for CSBND")}
            A=c(0:(m-j-1),(m-j+1):m,(v-(m-j)))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==4){
            j=(m-4)/8
            if(j<0) {return("Conditions are not satisfied for CSBND")}
            A=c(0:j,(j+2):(m-1),(m+1),(v-(j+1)))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==5){
            j=(m-5)/8
            if(j<0) {return("Conditions are not satisfied for CSBND")}
            A=c(0:(3*j+1),(3*j+3):(m),(v-(3*j+2)))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==6){
            j=(m-6)/8
            if(j<0) {return("Conditions are not satisfied for CSBND")}
            A=c(0:(5*j+3),(5*j+5):(m),(v-(5*j+4)))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==7){
            j=(m-7)/8
            if(j<1) {return("Conditions are not satisfied for CSBND")}
            A=c(0:(m-j-1),(m-j+1):(m-1),(m+1),(v-(m-j)))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          } 
        }  
        if(D==1 & sk2==2 & sk3==2){
          v= 2*i*k[1]+4*k[2]+4*k[3]+1  ;m=(v-1)/2
          if(m%%8==0){
            j=m/8
            if(j<1) {return("Conditions are not satisfied for CBNDs")}
            A=c(1:(j-1),(j+1):m,(v-j))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==1){
            j=(m-1)/8
            if(j<1) {return("Conditions are not satisfied for CBNDs")}
            A=c(1:(3*j),(3*j+2):(m-1),(m+1),(v-(3*j+1)))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==2){
            j=(m-2)/8
            if(j<1) {return("Conditions are not satisfied for CBNDs")}
            A=c(1:(5*j+1),(5*j+3):(m-1),(m+1),(v-(5*j+2)))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==3){
            j=(m-3)/8
            if(j<0) {return("Conditions are not satisfied for CBNDs")}
            A=c(1:(m-j-1),(m-j+1):m,(v-(m-j)))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==4){
            j=(m-4)/8
            if(j<0) {return("Conditions are not satisfied for CBNDs")}
            A=c(1:j,(j+2):(m-1),(m+1),(v-(j+1)))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==5){
            j=(m-5)/8
            if(j<0) {return("Conditions are not satisfied for CBNDs")}
            A=c(1:(3*j+1),(3*j+3):(m),(v-(3*j+2)))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==6){
            j=(m-6)/8
            if(j<0) {return("Conditions are not satisfied for CBNDs")}
            A=c(1:(5*j+3),(5*j+5):(m),(v-(5*j+4)))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==7){
            j=(m-7)/8
            if(j<1) {return("Conditions are not satisfied for CBNDs")}
            A=c(1:(m-j-1),(m-j+1):(m-1),(m+1),(v-(m-j)))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
        } 
        if(D==2 & sk2==1 & sk3==2){  
          v=2*i*k[1]+2*k[2]+4*k[3]-1 ; m=(v-1)/2
          if(m%%8==0){
            j=m/8
            if(j<1) {return("Conditions are not satisfied for CSBND")}
            A=c(0:(j-1),(j+1):m,(v-j))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==1){
            j=(m-1)/8
            if(j<1) {return("Conditions are not satisfied for CSBND")}
            A=c(0:(3*j),(3*j+2):(m-1),(m+1),(v-(3*j+1)))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==2){
            j=(m-2)/8
            if(j<1) {return("Conditions are not satisfied for CSBND")}
            A=c(0:(5*j+1),(5*j+3):(m-1),(m+1),(v-(5*j+2)))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==3){
            j=(m-3)/8
            if(j<0) {return("Conditions are not satisfied for CSBND")}
            A=c(0:(m-j-1),(m-j+1):m,(v-(m-j)))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==4){
            j=(m-4)/8
            if(j<0) {return("Conditions are not satisfied for CSBND")}
            A=c(0:j,(j+2):(m-1),(m+1),(v-(j+1)))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==5){
            j=(m-5)/8
            if(j<0) {return("Conditions are not satisfied for CSBND")}
            A=c(0:(3*j+1),(3*j+3):(m),(v-(3*j+2)))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==6){
            j=(m-6)/8
            if(j<0) {return("Conditions are not satisfied for CSBND")}
            A=c(0:(5*j+3),(5*j+5):(m),(v-(5*j+4)))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==7){
            j=(m-7)/8
            if(j<1) {return("Conditions are not satisfied for CSBND")}
            A=c(0:(m-j-1),(m-j+1):(m-1),(m+1),(v-(m-j)))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          } 
        }
        if(D==1 & sk2==1 & sk3==2){
          v= 2*i*k[1]+2*k[2]+4*k[3]+1  ;m=(v-1)/2
          if(m%%8==0){
            j=m/8
            if(j<1) {return("Conditions are not satisfied for CBNDs")}
            A=c(1:(j-1),(j+1):m,(v-j))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==1){
            j=(m-1)/8
            if(j<1) {return("Conditions are not satisfied for CBNDs")}
            A=c(1:(3*j),(3*j+2):(m-1),(m+1),(v-(3*j+1)))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==2){
            j=(m-2)/8
            if(j<1) {return("Conditions are not satisfied for CBNDs")}
            A=c(1:(5*j+1),(5*j+3):(m-1),(m+1),(v-(5*j+2)))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==3){
            j=(m-3)/8
            if(j<0) {return("Conditions are not satisfied for CBNDs")}
            A=c(1:(m-j-1),(m-j+1):m,(v-(m-j)))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==4){
            j=(m-4)/8
            if(j<0) {return("Conditions are not satisfied for CBNDs")}
            A=c(1:j,(j+2):(m-1),(m+1),(v-(j+1)))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==5){
            j=(m-5)/8
            if(j<0) {return("Conditions are not satisfied for CBNDs")}
            A=c(1:(3*j+1),(3*j+3):(m),(v-(3*j+2)))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==6){
            j=(m-6)/8
            if(j<0) {return("Conditions are not satisfied for CBNDs")}
            A=c(1:(5*j+3),(5*j+5):(m),(v-(5*j+4)))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==7){
            j=(m-7)/8
            if(j<1) {return("Conditions are not satisfied for CBNDs")}
            A=c(1:(m-j-1),(m-j+1):(m-1),(m+1),(v-(m-j)))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
        } 
        if(D==2 & sk2==2 & sk3==1){  
          v=2*i*k[1]+4*k[2]+2*k[3]-1; m=(v-1)/2
          if(m%%8==0){
            j=m/8
            if(j<1) {return("Conditions are not satisfied for CSBND")}
            A=c(0:(j-1),(j+1):m,(v-j))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==1){
            j=(m-1)/8
            if(j<1) {return("Conditions are not satisfied for CSBND")}
            A=c(0:(3*j),(3*j+2):(m-1),(m+1),(v-(3*j+1)))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==2){
            j=(m-2)/8
            if(j<1) {return("Conditions are not satisfied for CSBND")}
            A=c(0:(5*j+1),(5*j+3):(m-1),(m+1),(v-(5*j+2)))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==3){
            j=(m-3)/8
            if(j<0) {return("Conditions are not satisfied for CSBND")}
            A=c(0:(m-j-1),(m-j+1):m,(v-(m-j)))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==4){
            j=(m-4)/8
            if(j<0) {return("Conditions are not satisfied for CSBND")}
            A=c(0:j,(j+2):(m-1),(m+1),(v-(j+1)))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==5){
            j=(m-5)/8
            if(j<0) {return("Conditions are not satisfied for CSBND")}
            A=c(0:(3*j+1),(3*j+3):(m),(v-(3*j+2)))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==6){
            j=(m-6)/8
            if(j<0) {return("Conditions are not satisfied for CSBND")}
            A=c(0:(5*j+3),(5*j+5):(m),(v-(5*j+4)))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==7){
            j=(m-7)/8
            if(j<1) {return("Conditions are not satisfied for CSBND")}
            A=c(0:(m-j-1),(m-j+1):(m-1),(m+1),(v-(m-j)))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
        }  
        if(D==1 & sk2==2 & sk3==1){
          v= 2*i*k[1]+4*k[2]+2*k[3]+1  ;m=(v-1)/2
          if(m%%8==0){
            j=m/8
            if(j<1) {return("Conditions are not satisfied for CBNDs")}
            A=c(1:(j-1),(j+1):m,(v-j))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==1){
            j=(m-1)/8
            if(j<1) {return("Conditions are not satisfied for CBNDs")}
            A=c(1:(3*j),(3*j+2):(m-1),(m+1),(v-(3*j+1)))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==2){
            j=(m-2)/8
            if(j<1) {return("Conditions are not satisfied for CBNDs")}
            A=c(1:(5*j+1),(5*j+3):(m-1),(m+1),(v-(5*j+2)))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==3){
            j=(m-3)/8
            if(j<0) {return("Conditions are not satisfied for CBNDs")}
            A=c(1:(m-j-1),(m-j+1):m,(v-(m-j)))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==4){
            j=(m-4)/8
            if(j<0) {return("Conditions are not satisfied for CBNDs")}
            A=c(1:j,(j+2):(m-1),(m+1),(v-(j+1)))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==5){
            j=(m-5)/8
            if(j<0) {return("Conditions are not satisfied for CBNDs")}
            A=c(1:(3*j+1),(3*j+3):(m),(v-(3*j+2)))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==6){
            j=(m-6)/8
            if(j<0) {return("Conditions are not satisfied for CBNDs")}
            A=c(1:(5*j+3),(5*j+5):(m),(v-(5*j+4)))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==7){
            j=(m-7)/8
            if(j<1) {return("Conditions are not satisfied for CBNDs")}
            A=c(1:(m-j-1),(m-j+1):(m-1),(m+1),(v-(m-j)))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
        } 
        if(D==12 & sk2==1 & sk3==1){  
          v=2*i*k[1]+2*k[2]+2*k[3]-3 ; m=(v-1)/2
          if(m%%8==0){
            j=m/8
            if(j<1) {return("Conditions are not satisfied for CSBGND")}
            A=c(0:(j-1),(j+1),(j+2):(m-1),(m+1),(2*m-1-j),1)
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==1){
            j=(m-1)/8
            if(j<1) {return("Conditions are not satisfied for CSBGND")}
            A=c(0:(3*j),(3*j+2),(3*j+3):m,(2*m-3*j),1)
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==2){
            j=(m-2)/8
            if(j<1) {return("Conditions are not satisfied for CSBGND")}
            A=c(0:(5*j+1),(5*j+3),(5*j+4):m,(2*m-5*j-1),1)
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==3){
            j=(m-3)/8
            if(j<0) {return("Conditions are not satisfied for CSBGND")}
            A=c(0:(m-j),(m-j+2),(m-j+3):m,(m+j),1)
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==4){
            j=(m-4)/8
            if(j<0) {return("Conditions are not satisfied for CSBGND")}
            A=c(0:j,(j+2),(j+3):m,(2*m-j),1)
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==5){
            j=(m-5)/8
            if(j<0) {return("Conditions are not satisfied for CSBGND")}
            A=c(0:(3*j+2),(3*j+4),(3*j+5):(m-1),(m+1),(2*m-3*j-2),1)
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==6){
            j=(m-6)/8
            if(j<0) {return("Conditions are not satisfied for CSBGND")}
            A=c(0:(5*j+4),(5*j+6),(5*j+7):(m-1),(m+1),(2*m-5*j-4),1)
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==7){
            j=(m-7)/8
            if(j<1) {return("Conditions are not satisfied for CSBGND")}
            A=c(0:(m-j-1),(m-j+1),(m+j+2):m,(m+j+1),1)
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
        }  
        if(D==11 & sk2==1 & sk3==1){
          v= 2*i*k[1]+2*k[2]+2*k[3]+3; m=(v-1)/2
          if(m%%8==0){
            j=m/8
            if(j<1) {return("Conditions are not satisfied for CPBNDs")}
            A=c(2:j,(j+2),(j+3):(m-2),m,(m+2),(2*m-j))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==1){
            j=(m-1)/8
            if(j<1) {return("Conditions are not satisfied for CPBNDs")}
            A=c(2:(3*j-1),(3*j+1),(3*j+2):m,(2*m-(3*j-1)))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==2){
            j=(m-2)/8
            if(j<1) {return("Conditions are not satisfied for CPBNDs")}
            A=c(2:(5*j+2),(5*j+3):m,(2*m-(5*j)))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==3){
            j=(m-3)/8
            if(j<0) {return("Conditions are not satisfied for CPBNDs")}
            A=c(2:(m-j-1),(m-j+1),(m-j+2):(m+1),(m+2),(m+j+1))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==4){
            j=(m-4)/8
            if(j<0) {return("Conditions are not satisfied for CPBNDs")}
            A=c(2:(j+1),(j+3),(j+4):(m-2),(m+1),(m+2),(2*m-j-1))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          
          if(m%%8==5){
            j=(m-5)/8
            if(j<0) {return("Conditions are not satisfied for CPBNDs")}
            A=c(2:(3*j+1),(3*j+3),(3*j+4):(m-1),(m+1),(2*m-3*j-1))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==6){
            j=(m-6)/8
            if(j<0) {return("Conditions are not satisfied for CPBNDs")}
            A=c(2:(5*j+3),(5*j+5),(5*j+6):(m-1),(m+1),(2*m-5*j-3))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==7){
            j=(m-7)/8
            if(j<1) {return("Conditions are not satisfied for CPBNDs")}
            A=c(2:(m-j-2),(m-j),(m-j+1):m,(m+j+2))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          } 
          }
        if(D==12 & sk2==2 & sk3==2){  
          v=2*i*k[1]+4*k[2]+4*k[3]-3 ; m=(v-1)/2
          
          if(m%%8==0){
            j=m/8
            if(j<1) {return("Conditions are not satisfied for CSBGND")}
            A=c(0:(j-1),(j+1),(j+2):(m-1),(m+1),(2*m-1-j),1)
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==1){
            j=(m-1)/8
            if(j<1) {return("Conditions are not satisfied for CSBGND")}
            A=c(0:(3*j),(3*j+2),(3*j+3):m,(2*m-3*j),1)
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==2){
            j=(m-2)/8
            if(j<1) {return("Conditions are not satisfied for CSBGND")}
            A=c(0:(5*j+1),(5*j+3),(5*j+4):m,(2*m-5*j-1),1)
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          
          if(m%%8==3){
            j=(m-3)/8
            if(j<0) {return("Conditions are not satisfied for CSBGND")}
            A=c(0:(m-j),(m-j+2),(m-j+3):m,(m+j),1)
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==4){
            j=(m-4)/8
            if(j<0) {return("Conditions are not satisfied for CSBGND")}
            A=c(0:j,(j+2),(j+3):m,(2*m-j),1)
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          
          if(m%%8==5){
            j=(m-5)/8
            if(j<0) {return("Conditions are not satisfied for CSBGND")}
            A=c(0:(3*j+2),(3*j+4),(3*j+5):(m-1),(m+1),(2*m-3*j-2),1)
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==6){
            j=(m-6)/8
            if(j<0) {return("Conditions are not satisfied for CSBGND")}
            A=c(0:(5*j+4),(5*j+6),(5*j+7):(m-1),(m+1),(2*m-5*j-4),1)
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==7){
            j=(m-7)/8
            if(j<1) {return("Conditions are not satisfied for CSBGND")}
            A=c(0:(m-j-1),(m-j+1),(m+j+2):m,(m+j+1),1)
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          } 
        }  
        if(D==11 & sk2==2 & sk3==2){
          v= 2*i*k[1]+4*k[2]+4*k[3]+3  ;m=(v-1)/2
          if(m%%8==0){
            j=m/8
            if(j<1) {return("Conditions are not satisfied for CPBNDs")}
            A=c(2:j,(j+2),(j+3):(m-2),m,(m+2),(2*m-j))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==1){
            j=(m-1)/8
            if(j<1) {return("Conditions are not satisfied for CPBNDs")}
            A=c(2:(3*j-1),(3*j+1),(3*j+2):m,(2*m-(3*j-1)))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==2){
            j=(m-2)/8
            if(j<1) {return("Conditions are not satisfied for CPBNDs")}
            A=c(1:(5*j+1),(5*j+3):(m-1),(m+1),(v-(5*j+2)))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==3){
            j=(m-3)/8
            if(j<0) {return("Conditions are not satisfied for CPBNDs")}
            A=c(2:(m-j-1),(m-j+1),(m-j+2):(m+1),(m+2),(m+j+1))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==4){
            j=(m-4)/8
            if(j<0) {return("Conditions are not satisfied for CPBNDs")}
            A=c(2:(j+1),(j+3),(j+4):(m-2),(m+1),(m+2),(2*m-j-1))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==5){
            j=(m-5)/8
            if(j<0) {return("Conditions are not satisfied for CPBNDs")}
            A=c(2:(3*j+1),(3*j+3),(3*j+4):(m-1),(m+1),(2*m-3*j-1))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==6){
            j=(m-6)/8
            if(j<0) {return("Conditions are not satisfied for CPBNDs")}
            A=c(2:(5*j+3),(5*j+5),(5*j+6):(m-1),(m+1),(2*m-5*j-3))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==7){
            j=(m-7)/8
            if(j<1) {return("Conditions are not satisfied for CPBNDs")}
            A=c(2:(m-j-2),(m-j),(m-j+1):m,(m+j+2))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          } 
        } 
        if(D==12 & sk2==1 & sk3==2){  
          v=2*i*k[1]+2*k[2]+4*k[3]-3 ; m=(v-1)/2
          if(m%%8==0){
            j=m/8
            if(j<1) {return("Conditions are not satisfied for CSBGND")}
            A=c(0:(j-1),(j+1),(j+2):(m-1),(m+1),(2*m-1-j),1)
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==1){
            j=(m-1)/8
            if(j<1) {return("Conditions are not satisfied for CSBGND")}
            A=c(0:(3*j),(3*j+2),(3*j+3):m,(2*m-3*j),1)
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==2){
            j=(m-2)/8
            if(j<1) {return("Conditions are not satisfied for CSBGND")}
            A=c(0:(5*j+1),(5*j+3),(5*j+4):m,(2*m-5*j-1),1)
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==3){
            j=(m-3)/8
            if(j<0) {return("Conditions are not satisfied for CSBGND")}
            A=c(0:(m-j),(m-j+2),(m-j+3):m,(m+j),1)
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==4){
            j=(m-4)/8
            if(j<0) {return("Conditions are not satisfied for CSBGND")}
            A=c(0:j,(j+2),(j+3):m,(2*m-j),1)
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==5){
            j=(m-5)/8
            if(j<0) {return("Conditions are not satisfied for CSBGND")}
            A=c(0:(3*j+2),(3*j+4),(3*j+5):(m-1),(m+1),(2*m-3*j-2),1)
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==6){
            j=(m-6)/8
            if(j<0) {return("Conditions are not satisfied for CSBGND")}
            A=c(0:(5*j+4),(5*j+6),(5*j+7):(m-1),(m+1),(2*m-5*j-4),1)
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==7){
            j=(m-7)/8
            if(j<1) {return("Conditions are not satisfied for CSBGND")}
            A=c(0:(m-j-1),(m-j+1),(m+j+2):m,(m+j+1),1)
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          } 
        } 
        if(D==11 & sk2==1 & sk3==2){
          v= 2*i*k[1]+2*k[2]+4*k[3]+3  ;m=(v-1)/2
          if(m%%8==0){
            j=m/8
            if(j<1) {return("Conditions are not satisfied for CPBNDs")}
            A=c(2:j,(j+2),(j+3):(m-2),m,(m+2),(2*m-j))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==1){
            j=(m-1)/8
            if(j<1) {return("Conditions are not satisfied for CPBNDs")}
            A=c(2:(3*j-1),(3*j+1),(3*j+2):m,(2*m-(3*j-1)))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==2){
            j=(m-2)/8
            if(j<1) {return("Conditions are not satisfied for CPBNDs")}
            A=c(2:(5*j+2),(5*j+3):m,(2*m-(5*j)))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==3){
            j=(m-3)/8
            if(j<0) {return("Conditions are not satisfied for CPBNDs")}
            A=c(2:(m-j-1),(m-j+1),(m-j+2):(m+1),(m+2),(m+j+1))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==4){
            j=(m-4)/8
            if(j<0) {return("Conditions are not satisfied for CPBNDs")}
            A=c(2:(j+1),(j+3),(j+4):(m-2),(m+1),(m+2),(2*m-j-1))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==5){
            j=(m-5)/8
            if(j<0) {return("Conditions are not satisfied for CPBNDs")}
            A=c(2:(3*j+1),(3*j+3),(3*j+4):(m-1),(m+1),(2*m-3*j-1))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==6){
            j=(m-6)/8
            if(j<0) {return("Conditions are not satisfied for CPBNDs")}
            A=c(2:(5*j+3),(5*j+5),(5*j+6):(m-1),(m+1),(2*m-5*j-3))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==7){
            j=(m-7)/8
            if(j<1) {return("Conditions are not satisfied for CPBNDs")}
            A=c(2:(m-j-2),(m-j),(m-j+1):m,(m+j+2))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          } 
        } 
        if(D==12 & sk2==2 & sk3==1){  
          v=2*i*k[1]+4*k[2]+2*k[3]-3; m=(v-1)/2
          if(m%%8==0){
            j=m/8
            if(j<1) {return("Conditions are not satisfied for CSBGND")}
            A=c(0:(j-1),(j+1),(j+2):(m-1),(m+1),(2*m-1-j),1)
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==1){
            j=(m-1)/8
            if(j<1) {return("Conditions are not satisfied for CSBGND")}
            A=c(0:(3*j),(3*j+2),(3*j+3):m,(2*m-3*j),1)
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==2){
            j=(m-2)/8
            if(j<1) {return("Conditions are not satisfied for CSBGND")}
            A=c(0:(5*j+1),(5*j+3),(5*j+4):m,(2*m-5*j-1),1)
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==3){
            j=(m-3)/8
            if(j<0) {return("Conditions are not satisfied for CSBGND")}
            A=c(0:(m-j),(m-j+2),(m-j+3):m,(m+j),1)
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==4){
            j=(m-4)/8
            if(j<0) {return("Conditions are not satisfied for CSBGND")}
            A=c(0:j,(j+2),(j+3):m,(2*m-j),1)
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==5){
            j=(m-5)/8
            if(j<0) {return("Conditions are not satisfied for CSBGND")}
            A=c(0:(3*j+2),(3*j+4),(3*j+5):(m-1),(m+1),(2*m-3*j-2),1)
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==6){
            j=(m-6)/8
            if(j<0) {return("Conditions are not satisfied for CSBGND")}
            A=c(0:(5*j+4),(5*j+6),(5*j+7):(m-1),(m+1),(2*m-5*j-4),1)
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==7){
            j=(m-7)/8
            if(j<1) {return("Conditions are not satisfied for CSBGND")}
            A=c(0:(m-j-1),(m-j+1),(m+j+2):m,(m+j+1),1)
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          } 
        } 
        if(D==11 & sk2==2 & sk3==1){
          v= 2*i*k[1]+4*k[2]+2*k[3]+3  ;m=(v-1)/2
          if(m%%8==0){
            j=m/8
            if(j<1) {return("Conditions are not satisfied for CPBNDs")}
            A=c(2:j,(j+2),(j+3):(m-2),m,(m+2),(2*m-j))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==1){
            j=(m-1)/8
            if(j<1) {return("Conditions are not satisfied for CPBNDs")}
            A=c(2:(3*j-1),(3*j+1),(3*j+2):m,(2*m-(3*j-1)))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==2){
            j=(m-2)/8
            if(j<1) {return("Conditions are not satisfied for CPBNDs")}
            A=c(2:(5*j+2),(5*j+3):m,(2*m-(5*j)))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==3){
            j=(m-3)/8
            if(j<0) {return("Conditions are not satisfied for CPBNDs")}
            A=c(2:(m-j-1),(m-j+1),(m-j+2):(m+1),(m+2),(m+j+1))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==4){
            j=(m-4)/8
            if(j<0) {return("Conditions are not satisfied for CPBNDs")}
            A=c(2:(j+1),(j+3),(j+4):(m-2),(m+1),(m+2),(2*m-j-1))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==5){
            j=(m-5)/8
            if(j<0) {return("Conditions are not satisfied for CPBNDs")}
            A=c(2:(3*j+1),(3*j+3),(3*j+4):(m-1),(m+1),(2*m-3*j-1))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==6){
            j=(m-6)/8
            if(j<0) {return("Conditions are not satisfied for CPBNDs")}
            A=c(2:(5*j+3),(5*j+5),(5*j+6):(m-1),(m+1),(2*m-5*j-3))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%8==7){
            j=(m-7)/8
            if(j<1) {return("Conditions are not satisfied for CPBNDs")}
            A=c(2:(m-j-2),(m-j),(m-j+1):m,(m+j+2))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          } 
        }
        if(D==3 & sk2==1 & sk3==1){  
          v=2*i*k[1]+2*k[2]+2*k[3]+2 ; m=(v-2)/2
          if(m%%4==0){
            A<-1:m
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%4==3){
            A<-c(1:((3*m-1)/4),((3*m+7)/4),((3*m+11)/4):m, (5*(m+1)/4))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%4==1 | m%%4==2){return("The minimal CPNDs in which v/2 unordered pair cannot be constructed for V=2ik1+2k2+2k3+2 and k1,k2 k3 block sizes")}
        } 
        if(D==4 & sk2==1 & sk3==1){
          v= 2*i*k[1]+2*k[2]+2*k[3]+4 
          m<-(v-2)/2
          if(m%%4==0){
            A=c(1:(m/2),((m+4)/2),((m+6)/2):(m-1),((3*m+2)/2))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%4==1){
            A=c(1:((3*m+1)/4),((3*m+9)/4),((3*m+13)/4):(m-1),((5*m+3)/4))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%4==2){
            A=c(2:(m-2),m,(2*m+1))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%4==3){
            A=c(1:((m+1)/4),((m+9)/4),((m+13)/4):(m-2), m, ((7*m+3)/4))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
        }
        if(D==3 & sk2==2 & sk3==2){  
          v=2*i*k[1]+4*k[2]+4*k[3]+2 ; m=(v-2)/2
          if(m%%4==0){
            A<-1:m
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%4==3){
            A<-c(1:((3*m-1)/4),((3*m+7)/4),((3*m+11)/4):m, (5*(m+1)/4))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%4==1 | m%%4==2){return("The minimal CPNDs in which v/2 unordered pair cannot be constructed for V=2ik1+4k2+4k3+2 and k1,k2 k3 block sizes")}
        }  
        if(D==4 & sk2==2 & sk3==2){
          v= 2*i*k[1]+4*k[2]+4*k[3]+4  ;m<-(v-2)/2
          if(m%%4==0){
            A=c(1:(m/2),((m+4)/2),((m+6)/2):(m-1),((3*m+2)/2))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%4==1){
            A=c(1:((3*m+1)/4),((3*m+9)/4),((3*m+13)/4):(m-1),((5*m+3)/4))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%4==2){
            A=c(2:(m-2),m,(2*m+1))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%4==3){
            A=c(1:((m+1)/4),((m+9)/4),((m+13)/4):(m-2), m, ((7*m+3)/4))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
        } 
        if(D==3 & sk2==1 & sk3==2){  
          v=2*i*k[1]+2*k[2]+4*k[3]+2 ; m=(v-2)/2
          if(m%%4==0){
            A<-1:m
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%4==3){
            A<-c(1:((3*m-1)/4),((3*m+7)/4),((3*m+11)/4):m, (5*(m+1)/4))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%4==1 | m%%4==2){return("The minimal CPNDs in which v/2 unordered pair cannot be constructed for V=2ik1+2k2+4k3+2 and k1,k2 k3 block sizes")}
        }  
        if(D==4 & sk2==1 & sk3==2){
          v= 2*i*k[1]+2*k[2]+4*k[3]+4  ;m<-(v-2)/2
          if(m%%4==0){
            A=c(1:(m/2),((m+4)/2),((m+6)/2):(m-1),((3*m+2)/2))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%4==1){
            A=c(1:((3*m+1)/4),((3*m+9)/4),((3*m+13)/4):(m-1),((5*m+3)/4))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%4==2){
            A=c(2:(m-2),m,(2*m+1))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%4==3){
            A=c(1:((m+1)/4),((m+9)/4),((m+13)/4):(m-2), m, ((7*m+3)/4))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
        } 
        if(D==3 & sk2==2 & sk3==1){  
          v=2*i*k[1]+4*k[2]+2*k[3]+2 ; m=(v-2)/2
          if(m%%4==0){
            A<-1:m
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%4==3){
            A<-c(1:((3*m-1)/4),((3*m+7)/4),((3*m+11)/4):m, (5*(m+1)/4))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%4==1 | m%%4==2){return("The minimal CPNDs in which v/2 unordered pair cannot be constructed for V=2ik1+4k2+2k3+2 and k1,k2 k3 block sizes")}
        }  
        if(D==4 & sk2==2 & sk3==1){
          v= 2*i*k[1]+4*k[2]+2*k[3]+4  ;m<-(v-2)/2
          if(m%%4==0){
            A=c(1:(m/2),((m+4)/2),((m+6)/2):(m-1),((3*m+2)/2))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%4==1){
            A=c(1:((3*m+1)/4),((3*m+9)/4),((3*m+13)/4):(m-1),((5*m+3)/4))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%4==2){
            A=c(2:(m-2),m,(2*m+1))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%4==3){
            A=c(1:((m+1)/4),((m+9)/4),((m+13)/4):(m-2), m, ((7*m+3)/4))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
        } 
        if(D==5 & sk2==1 & sk3==1){  
          v=2*i*k[1]+2*k[2]+2*k[3] ; m<-(v-2)/2
          if(m%%4==2){
            A=1:(m+1)
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%4==3){
            A<-c(1:((m-3)/4),((m+5)/4),((m+9)/4):m,(m+1),(7*(m+1)/4))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%4==0 | m%%4==1){return("The minimal CWBNDs in which v/2 unordered pairs apears twice cannot be constructed for v=2ik1+2k2+2k3 and and k1,k2 k3 block sizes")}
        }  
        if(D==6 & sk2==1 & sk3==1){
          v= 2*i*k[1]+2*k[2]+2*k[3]-2; m<-(v-2)/2
          if(m%%4==0){
            A=c(1:(m-1),(m-1),(m+1),(m+2))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%4==1){
            A=c(1:((m-5)/4),((m+3)/4),((m+7)/4):m,(m+1),m,((7*m+9)/4))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%4==2){
            A=c(1:((m-2)/2),((m+2)/2),((m+4)/2):m,(m+1),m,((3*m+4)/2))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%4==3){
            A=c(1:((3*m-5)/4),((3*m+3)/4),((3*m+7)/4):(m-1), m,(m+1),(m-1), ((5*m+9)/4))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
        }   
        if(D==5 & sk2==2 & sk3==2){  
          v=2*i*k[1]+4*k[2]+4*k[3] ; m=(v-2)/2
          if(m%%4==2){
            A=1:(m+1)
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%4==3){
            A<-c(1:((m-3)/4),((m+5)/4),((m+9)/4):m,(m+1),(7*(m+1)/4))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%4==0 | m%%4==1){return("The minimal CWBNDs in which v/2 unordered pairs apears twice cannot be constructed for v=2ik1+4k2+4k3 and and k1,k2 k3 block sizes")}
        }  
        if(D==6 & sk2==2 & sk3==2){
          v= 2*i*k[1]+4*k[2]+4*k[3]-2  ;m<-(v-2)/2
          if(m%%4==0){
            A=c(1:(m-1),(m-1),(m+1),(m+2))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%4==1){
            A=c(1:((m-5)/4),((m+3)/4),((m+7)/4):m,(m+1),m,((7*m+9)/4))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%4==2){
            A=c(1:((m-2)/2),((m+2)/2),((m+4)/2):m,(m+1),m,((3*m+4)/2))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%4==3){
            A=c(1:((3*m-5)/4),((3*m+3)/4),((3*m+7)/4):(m-1), m,(m+1),(m-1), ((5*m+9)/4))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
        } 
        if(D==5 & sk2==1 & sk3==2){  
          v=2*i*k[1]+2*k[2]+4*k[3] ; m=(v-2)/2
          if(m%%4==2){
            A=1:(m+1)
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%4==3){
            A<-c(1:((m-3)/4),((m+5)/4),((m+9)/4):m,(m+1),(7*(m+1)/4))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%4==0 | m%%4==1){return("The minimal CWBNDs in which v/2 unordered pairs apears twice cannot be constructed for v=2ik1+2k2+4k3 and and k1,k2 k3 block sizes")}
        }  
        if(D==6 & sk2==1 & sk3==2){
          v= 2*i*k[1]+2*k[2]+4*k[3]-2  ;m<-(v-2)/2
          if(m%%4==0){
            A=c(1:(m-1),(m-1),(m+1),(m+2))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%4==1){
            A=c(1:((m-5)/4),((m+3)/4),((m+7)/4):m,(m+1),m,((7*m+9)/4))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%4==2){
            A=c(1:((m-2)/2),((m+2)/2),((m+4)/2):m,(m+1),m,((3*m+4)/2))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%4==3){
            A=c(1:((3*m-5)/4),((3*m+3)/4),((3*m+7)/4):(m-1), m,(m+1),(m-1), ((5*m+9)/4))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
        } 
        if(D==5 & sk2==2 & sk3==1){  
          v=2*i*k[1]+4*k[2]+2*k[3] ; m=(v-2)/2
          if(m%%4==2){
            A=1:(m+1)
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%4==3){
            A<-c(1:((m-3)/4),((m+5)/4),((m+9)/4):m,(m+1),(7*(m+1)/4))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%4==0 | m%%4==1){return("The minimal CWBNDs in which v/2 unordered pairs apears twice cannot be constructed for v=2ik1+4k2+2k3 and and k1,k2 k3 block sizes")}
        }  
        if(D==6 & sk2==2 & sk3==1){
          v= 2*i*k[1]+4*k[2]+2*k[3]-2  ;m<-(v-2)/2
          if(m%%4==0){
            A=c(1:(m-1),(m-1),(m+1),(m+2))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%4==1){
            A=c(1:((m-5)/4),((m+3)/4),((m+7)/4):m,(m+1),m,((7*m+9)/4))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%4==2){
            A=c(1:((m-2)/2),((m+2)/2),((m+4)/2):m,(m+1),m,((3*m+4)/2))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%4==3){
            A=c(1:((3*m-5)/4),((3*m+3)/4),((3*m+7)/4):(m-1), m,(m+1),(m-1), ((5*m+9)/4))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
        }
        if(D==7 & sk2==1 & sk3==1){ 
          v=2*i*k[1]+2*k[2]+2*k[3] ; m=(v-2)/2
          if(m%%4==0){
            A<-0:m
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%4==3){
            A<-c(0:((3*m-1)/4),((3*m+7)/4),((3*m+11)/4):m, (5*(m+1)/4))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%4==1 | m%%4==2){return("The minimal CSPBNDs in which v/2 unordered pair cannot be constructed for V=2ik1+2k2+2k3+2 and k1,k2 k3 block sizes")}
        } 
        if(D==8 & sk2==1 & sk3==1){
          v= 2*i*k[1]+2*k[2]+2*k[3]+2
          m<-(v-2)/2
          if(m%%4==0){
            A=c(0:(m/2),((m+4)/2),((m+6)/2):(m-1),((3*m+2)/2))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%4==1){
            A=c(0:((3*m+1)/4),((3*m+9)/4),((3*m+13)/4):(m-1),((5*m+3)/4))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%4==2){
            A=c(0,2:(m-2),m,(2*m+1))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%4==3){
            A=c(0:((m+1)/4),((m+9)/4),((m+13)/4):(m-2), m, ((7*m+3)/4))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
        } 
        if(D==7 & sk2==2 & sk3==2){ 
          v=2*i*k[1]+4*k[2]+4*k[3] ; m=(v-2)/2
          if(m%%4==0){
            A<-0:m
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%4==3){
            A<-c(0:((3*m-1)/4),((3*m+7)/4),((3*m+11)/4):m, (5*(m+1)/4))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%4==1 | m%%4==2){return("The minimal CSPBNDs in which v/2 unordered pair cannot be constructed for V=2ik1+4k2+4k3+2 and k1,k2 k3 block sizes")}
        }
        if(D==8 & sk2==2 & sk3==2){
          v= 2*i*k[1]+4*k[2]+4*k[3]+2 ;m<-(v-2)/2
          if(m%%4==0){
            A=c(0:(m/2),((m+4)/2),((m+6)/2):(m-1),((3*m+2)/2))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%4==1){
            A=c(0:((3*m+1)/4),((3*m+9)/4),((3*m+13)/4):(m-1),((5*m+3)/4))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%4==2){
            A=c(0,2:(m-2),m,(2*m+1))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%4==3){
            A=c(0:((m+1)/4),((m+9)/4),((m+13)/4):(m-2), m, ((7*m+3)/4))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
        } 
        if(D==7 & sk2==1 & sk3==2){ 
          v=2*i*k[1]+2*k[2]+4*k[3] ; m=(v-2)/2
          if(m%%4==0){
            A<-0:m
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%4==3){
            A<-c(0:((3*m-1)/4),((3*m+7)/4),((3*m+11)/4):m, (5*(m+1)/4))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%4==1 | m%%4==2){return("The minimal CSPBNDs in which v/2 unordered pair cannot be constructed for V=2ik1+2k2+4k3+2 and k1,k2 k3 block sizes")}
        } 
        if(D==8 & sk2==1 & sk3==2){
          v= 2*i*k[1]+2*k[2]+4*k[3]+2 ;m<-(v-2)/2
          if(m%%4==0){
            A=c(0:(m/2),((m+4)/2),((m+6)/2):(m-1),((3*m+2)/2))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%4==1){
            A=c(0:((3*m+1)/4),((3*m+9)/4),((3*m+13)/4):(m-1),((5*m+3)/4))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%4==2){
            A=c(0,2:(m-2),m,(2*m+1))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%4==3){
            A=c(0:((m+1)/4),((m+9)/4),((m+13)/4):(m-2), m, ((7*m+3)/4))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
        } 
        if(D==7 & sk2==2 & sk3==1){ 
          v=2*i*k[1]+4*k[2]+2*k[3] ; m=(v-2)/2
          if(m%%4==0){
            A<-0:m
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%4==3){
            A<-c(0:((3*m-1)/4),((3*m+7)/4),((3*m+11)/4):m, (5*(m+1)/4))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%4==1 | m%%4==2){return("The minimal CSPBNDs in which v/2 unordered pair cannot be constructed for V=2ik1+4k2+2k3+2 and k1,k2 k3 block sizes")}
        } 
        if(D==8 & sk2==2 & sk3==1){
          v= 2*i*k[1]+4*k[2]+2*k[3]+2 ;m<-(v-2)/2
          if(m%%4==0){
            A=c(0:(m/2),((m+4)/2),((m+6)/2):(m-1),((3*m+2)/2))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%4==1){
            A=c(0:((3*m+1)/4),((3*m+9)/4),((3*m+13)/4):(m-1),((5*m+3)/4))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%4==2){
            A=c(0,2:(m-2),m,(2*m+1))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%4==3){
            A=c(0:((m+1)/4),((m+9)/4),((m+13)/4):(m-2), m, ((7*m+3)/4))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
        }
        if(D==9 & sk2==1 & sk3==1){  
          v=2*i*k[1]+2*k[2]+2*k[3]-2 ; m<-(v-2)/2
          if(m%%4==2){
            A=0:(m+1)
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%4==3){
            A<-c(0:((m-3)/4),((m+5)/4),((m+9)/4):m,(m+1),(7*(m+1)/4))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%4==0 | m%%4==1){return("The minimal CSBGNDs in which v/2 unordered pairs apears twice cannot be constructed for v=2ik1+2k2+2k3 and and k1,k2 k3 block sizes")}
        }
        if(D==10 & sk2==1 & sk3==1){
          v= 2*i*k[1]+2*k[2]+2*k[3]-4; m<-(v-2)/2
          if(m%%4==0){
            A=c(0:(m-1),(m-1),(m+1),(m+2))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%4==1){
            A=c(0:((m-5)/4),((m+3)/4),((m+7)/4):m,(m+1),m,((7*m+9)/4))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%4==2){
            A=c(0:((m-2)/2),((m+2)/2),((m+4)/2):m,(m+1),m,((3*m+4)/2))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%4==3){
            A=c(0:((3*m-5)/4),((3*m+3)/4),((3*m+7)/4):(m-1), m,(m+1),(m-1), ((5*m+9)/4))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }  
        }  
        if(D==9 & sk2==2 & sk3==2){  
          v=2*i*k[1]+4*k[2]+4*k[3]-2 ; m=(v-2)/2
          if(m%%4==2){
            A=0:(m+1)
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%4==3){
            A<-c(0:((m-3)/4),((m+5)/4),((m+9)/4):m,(m+1),(7*(m+1)/4))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%4==0 | m%%4==1){return("The minimal CSBGNDs in which v/2 unordered pairs apears twice cannot be constructed for v=2ik1+4k2+4k3 and and k1,k2 k3 block sizes")}
        }
        if(D==10 & sk2==2 & sk3==2){
          v= 2*i*k[1]+4*k[2]+4*k[3]-4  ;m<-(v-2)/2
          if(m%%4==0){
            A=c(0:(m-1),(m-1),(m+1),(m+2))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%4==1){
            A=c(0:((m-5)/4),((m+3)/4),((m+7)/4):m,(m+1),m,((7*m+9)/4))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%4==2){
            A=c(0:((m-2)/2),((m+2)/2),((m+4)/2):m,(m+1),m,((3*m+4)/2))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%4==3){
            A=c(0:((3*m-5)/4),((3*m+3)/4),((3*m+7)/4):(m-1), m,(m+1),(m-1), ((5*m+9)/4))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
        }
        if(D==9 & sk2==1 & sk3==2){  
          v=2*i*k[1]+2*k[2]+4*k[3]-2 ; m=(v-2)/2
          if(m%%4==2){
            A=0:(m+1)
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%4==3){
            A<-c(0:((m-3)/4),((m+5)/4),((m+9)/4):m,(m+1),(7*(m+1)/4))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%4==0 | m%%4==1){return("The minimal CSBGNDs in which v/2 unordered pairs apears twice cannot be constructed for v=2ik1+2k2+4k3 and and k1,k2 k3 block sizes")}
        } 
        if(D==10 & sk2==1 & sk3==2){
          v= 2*i*k[1]+2*k[2]+4*k[3]-4  ;m<-(v-2)/2
          if(m%%4==0){
            A=c(0:(m-1),(m-1),(m+1),(m+2))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%4==1){
            A=c(0:((m-5)/4),((m+3)/4),((m+7)/4):m,(m+1),m,((7*m+9)/4))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%4==2){
            A=c(0:((m-2)/2),((m+2)/2),((m+4)/2):m,(m+1),m,((3*m+4)/2))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%4==3){
            A=c(0:((3*m-5)/4),((3*m+3)/4),((3*m+7)/4):(m-1), m,(m+1),(m-1), ((5*m+9)/4))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
        }
        if(D==9 & sk2==2 & sk3==1){  
          v=2*i*k[1]+4*k[2]+2*k[3]-2 ; m=(v-2)/2
          if(m%%4==2){
            A=0:(m+1)
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%4==3){
            A<-c(0:((m-3)/4),((m+5)/4),((m+9)/4):m,(m+1),(7*(m+1)/4))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%4==0 | m%%4==1){return("The minimal CSBGNDs in which v/2 unordered pairs apears twice cannot be constructed for v=2ik1+4k2+2k3 and and k1,k2 k3 block sizes")}
        }
        if(D==10 & sk2==2 & sk3==1){
          v= 2*i*k[1]+4*k[2]+2*k[3]-4  ;m<-(v-2)/2
          if(m%%4==0){
            A=c(0:(m-1),(m-1),(m+1),(m+2))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%4==1){
            A=c(0:((m-5)/4),((m+3)/4),((m+7)/4):m,(m+1),m,((7*m+9)/4))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%4==2){
            A=c(0:((m-2)/2),((m+2)/2),((m+4)/2):m,(m+1),m,((3*m+4)/2))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
          if(m%%4==3){
            A=c(0:((3*m-5)/4),((3*m+3)/4),((3*m+7)/4):(m-1), m,(m+1),(m-1), ((5*m+9)/4))
            A1<-grouping3(A,k,v,i,sk2,sk3)
            A2<-c(v,k);names(A2)<-c("V","K1","K2","K3")
            x<-list(S=A1$B1,G=A1$B4,R=A2,A=A)
          }
        } 
        new("stat_test", x)  
      }
      ###################################################################
      # Generation of design using sets of cyclical shifts
            design_CBND<-function(H){
        
        setClass( "CBND_design", representation("list"))
        setMethod("show", "CBND_design", function(object) {
          row <- paste(rep("=", 51), collapse = "")
          cat(row, "\n")
          cat("Following is minimal CBND for", "v=" ,object$R[1], "and","k=",object$R[2], "\n")
          row <- paste(rep("=", 51), collapse = "")
          cat(row, "\n")
          for(i in 1:length(ss)){
            W<-ss[[i]]
            nr<-dim(W)[1]
            for(j in 1:nr){
              print(object$Design[[i]][[j]])
              cat("\n\n")
            }}
        })  
        v<-H$R[1]
        k<-H$R[2]
        ss<-H$S  
        treat<-(1:v)-1
        fn<-(1:v)
        G<-list()
        for(j in 1:length(ss)){ 
          W<-ss[[j]]
          nr<-dim(W)[1]
          nc<-dim(W)[2]
          D<-list()
          for(i in 1:nr){
            dd<-c()
            d1<-matrix(treat,(nc+1),v,byrow = T)
            ss1<-cumsum(c(0,W[i,]))
            dd2<-d1+ss1
            dd<-rbind(dd,dd2)
            rr<-dd[which(dd>=v)]%%v
            dd[which(dd>=v)]<-rr
            colnames(dd)<-paste("B",fn, sep="")
            rownames(dd)<-rep("",(nc+1))
            fn<-fn+v
            D[[i]]<-dd
          }
          G[[j]]<-D
        }
        x<-list(Design=G,R=H$R)
        new("CBND_design", x)
      }
      Designs<-function(v,k1,k2,k3,i){
        if(v%%2!=0){
          if((v-1-(2*k2)-(2*k3))%%(2*k1)==0 & (i=(v-1-(2*k2)-(2*k3))/(2*k1))& (i%%1==0) &(i>0)){
            cat("MCBNDs is possible equale for both k2 and k3 ","v=",v,"k1=",k1,"k2=",k2,"k3=",k3, "i=", i )
            print(H<-CBND_3diffsize(k=c(k1,k2,k3),i=i,D=1,sk2=1,sk3=1))
          }else
            if((v-1-(4*k2)-(2*k3))%%(2*k1)==0 & (i=(v-1-(4*k2)-(2*k3))/(2*k1))& (i%%1==0) &(i>0)){
              cat("MCBNDs is possible with two sets for  k2 and one set for k3 ","v=",v,"k1=",k1,"k2=",k2,"k3=",k3, "i=", i )
              print(H<-CBND_3diffsize(k=c(k1,k2,k3),i=i,D=1,sk2=2,sk3=1)) 
            }else
              if((v-1-(2*k2)-(4*k3))%%(2*k1)==0 & (i=(v-1-(2*k2)-(4*k3))/(2*k1))& (i%%1==0) &(i>0)){
                
                cat("MCBNDs is possible with one set for  k2 and  two sets for k3 ","v=",v,"k1=",k1,"k2=",k2,"k3=",k3, "i=", i )
                print(H<-CBND_3diffsize(k=c(k1,k2,k3),i=i,D=1,sk2=1,sk3=2))
              }else
                if((v-1-(4*k2)-(4*k3))%%(2*k1)==0 & (i=(v-1-(4*k2)-(2*k3))/(2*k1))& (i%%1==0) &(i>0)){
                  cat("MCBNDs is possible with two sets for k2 and k3 ","v=",v,"k1=",k1,"k2=",k2,"k3=",k3,"k3=",k3,"i=",i )
                  print(H<-CBND_3diffsize(k=c(k1,k2,k3),i=i,D=1,sk2=2,sk3=2))
                }
                  if(((v+1-(2*k2)-(2*k3))%%(2*k1)==0)& (i=(v+1-(2*k2)-(2*k3))/(2*k1))& (i%%1==0) &(i>0)){
                    cat("MCSBNDs is possible  one set for k2 and k3","v=",v,"k1=",k1,"k2=",k2,"k3=",k3,"i=",i )
                    print(H<-CBND_3diffsize(k=c(k1,k2,k3),i=i,D=2,sk2=1,sk3=1))
                  }else
                    if((v+1-(4*k2)-(2*k3))%%(2*k1)==0 & (i=(v+1-(4*k2)-(2*k3))/(2*k1))& (i%%1==0) &(i>0)){
                      cat("MCSBNDs is possible with two set of k2 and one set for k3","v=",v,"k1=",k1,"k2=",k2,"k3=",k3,"i=",i )
                      print(H<-CBND_3diffsize(k=c(k1,k2,k3),i=i,D=2,sk2=2,sk3=1))
                    }else
                      if((v+1-(2*k2)-(4*k3))%%(2*k1)==0 & (i=(v+1-(2*k2)-(4*k3))/(2*k1))& (i%%1==0) &(i>0)){
                        cat("MCSBNDs is possible with one set of k2 and two set for k3","v=",v,"k1=",k1,"k2=",k2,"k3=",k3,"i=",i )
                        print(H<-CBND_3diffsize(k=c(k1,k2,k3),i=i,D=2,sk2=1,sk3=2))
                      }else
                        if((v+1-(4*k2)-(4*k3))%%(2*k1)==0 & (i=(v+1-(4*k2)-(4*k3))/(2*k1))& (i%%1==0) &(i>0)){
                          cat("MCSBNDs is possible with two set of k2 and k3","v=",v,"k1=",k1,"k2=",k2,"k3=",k3,"i=",i )
                          print(H<-CBND_3diffsize(k=c(k1,k2,k3),i=i,D=2,sk2=2,sk3=2))
                        }
                          if((v-3-(2*k2)-(2*k3))%%(2*k1)==0 & (i=(v-3-(2*k2)-(2*k3))/(2*k1))& (i%%1==0) &(i>0)){
                            cat("MCPBND-Odd is possible equale for both k2 and k3 ","v=",v,"k1=",k1,"k2=",k2,"k3=",k3, "i=", i )
                            print(H<-CBND_3diffsize(k=c(k1,k2,k3),i=i,D=11,sk2=1,sk3=1))
                          }else
                            if((v-3-(4*k2)-(2*k3))%%(2*k1)==0 & (i=(v-3-(4*k2)-(2*k3))/(2*k1))& (i%%1==0) &(i>0)){
                              cat("MCPBND-Odd is possible with two sets for  k2 and one set for k3 ","v=",v,"k1=",k1,"k2=",k2,"k3=",k3, "i=", i )
                              print(H<-CBND_3diffsize(k=c(k1,k2,k3),i=i,D=11,sk2=2,sk3=1)) 
                            }else
                              if((v-3-(2*k2)-(4*k3))%%(2*k1)==0 & (i=(v-3-(2*k2)-(4*k3))/(2*k1))& (i%%1==0) &(i>0)){
                                cat("MCPBND-Odd is possible with one set for  k2 and  two sets for k3 ","v=",v,"k1=",k1,"k2=",k2,"k3=",k3, "i=", i )
                                print(H<-CBND_3diffsize(k=c(k1,k2,k3),i=i,D=11,sk2=1,sk3=2))
                              }else
                                if((v-3-(4*k2)-(4*k3))%%(2*k1)==0 & (i=(v-3-(4*k2)-(2*k3))/(2*k1))& (i%%1==0) &(i>0)){
                                  cat("MCPBND-Odd is possible with two sets for k2 and k3 ","v=",v,"k1=",k1,"k2=",k2,"k3=",k3,"k3=",k3,"i=",i )
                                  print(H<-CBND_3diffsize(k=c(k1,k2,k3),i=i,D=11,sk2=2,sk3=2))
                                }
                                  if(((v+3-(2*k2)-(2*k3))%%(2*k1)==0)& (i=(v+3-(2*k2)-(2*k3))/(2*k1))& (i%%1==0) &(i>0)){
                                    cat("MCSBGND-Odd is possible  one set for k2 and k3","v=",v,"k1=",k1,"k2=",k2,"k3=",k3,"i=",i )
                                    print(H<-CBND_3diffsize(k=c(k1,k2,k3),i=i,D=12,sk2=1,sk3=1))
                                  }else
                                    if((v+3-(4*k2)-(2*k3))%%(2*k1)==0 & (i=(v+3-(4*k2)-(2*k3))/(2*k1))& (i%%1==0) &(i>0)){
                                      cat("MCSBGND-Odd is possible with two set of k2 and one set for k3","v=",v,"k1=",k1,"k2=",k2,"k3=",k3,"i=",i )
                                      print(H<-CBND_3diffsize(k=c(k1,k2,k3),i=i,D=12,sk2=2,sk3=1))
                                    }else
                                      if((v+3-(2*k2)-(4*k3))%%(2*k1)==0 & (i=(v+3-(2*k2)-(4*k3))/(2*k1))& (i%%1==0) &(i>0)){
                                        cat("MCSBGND-Odd is possible with one set of k2 and two set for k3","v=",v,"k1=",k1,"k2=",k2,"k3=",k3,"i=",i )
                                        print(H<-CBND_3diffsize(k=c(k1,k2,k3),i=i,D=12,sk2=1,sk3=2))
                                      }else
                                        if((v+3-(4*k2)-(4*k3))%%(2*k1)==0 & (i=(v+3-(4*k2)-(4*k3))/(2*k1))& (i%%1==0) &(i>0)){
                                          cat("MSBGND-Odd is possible with two set of k2 and k3","v=",v,"k1=",k1,"k2=",k2,"k3=",k3,"i=",i )
                                          print(H<-CBND_3diffsize(k=c(k1,k2,k3),i=i,D=12,sk2=2,sk3=2))
                                        }
          if((i%%1<0) |(i%%1>0))cat("MCBNDs,MCSBNDs, MCPBND-Odd and MSBGND-Odd are not possible","v=",v,"k1=",k1,"k2=",k2,"k3=",k3,"i=",i)
        }else
          if(v%%2==0){
            if((v-(2*k2)-(2*k3)-2)%%(2*k1)==0 &(i=(v-(2*k2)-(2*k3)-2)/(2*k1)) & (i%%1==0) &(i>0) & (((v-2)/2)%%4==0 |((v-2)/2)%%4==3)){
              cat("MCPBNDs-I is possible with one set of k2 and k3","v=",v,"k1=",k1,"k2=",k2,"k3=",k3,"i=",i )
             print (H<-CBND_3diffsize(k=c(k1,k2,k3),i=i,D=3,sk2=1,sk3=1))
            }else
              if((v-(4*k2)-(2*k3)-2)%%(2*k1)==0 &(i=(v-(4*k2)-(2*k3)-2)/(2*k1))&  (i%%1==0) &(i>0) & (((v-2)/2)%%4==0 |((v-2)/2)%%4==3)){
                cat("MCPBNDs-I is possible with two sets of k2 and one set of k3","v=",v,"k1=",k1,"k2=",k2,"k3=",k3,"i=",i )
                print(H<-CBND_3diffsize(k=c(k1,k2,k3),i=i,D=3,sk2=2,sk3=1))
              }else
                if((v-(2*k2)-(4*k3)-2)%%(2*k1)==0&(i=(v-(2*k2)-(4*k3)-2)/(2*k1))&  (i%%1==0) &(i>0) & (((v-2)/2)%%4==0 |((v-2)/2)%%4==3)){
                  cat("MCPBNDs-I is possible with one set of k2 and two sets of k3","v=",v,"k1=",k1,"k2=",k2,"k3=",k3,"i=",i )
                  print(H<-CBND_3diffsize(k=c(k1,k2,k3),i=i,D=3,sk2=1,sk3=2))
                }else
                  if((v-(4*k2)-(4*k3)-2)%%(2*k1)==0& (i=(v-(4*k2)-(4*k3)-2)/(2*k1)) &(i%%1==0) &(i>0) & (((v-2)/2)%%4==0 |((v-2)/2)%%4==3)){
                    cat("MCPBNDs-I is possible with two sets of k2 and k3","v=",v,"k1=",k1,"k2=",k2,"k3=",k3,"i=",i )
                    print(H<-CBND_3diffsize(k=c(k1,k2,k3),i=i,D=3,sk2=2,sk3=2))
                  }
                    if((v-(2*k2)-(2*k3)-2)%%(2*k1)==0 & (i=(v-(2*k2)-(2*k3)-2)/(2*k1)) & (i%%1==0) &(i>0) ){
                      cat("MCSPBNDS-II is possible with one set of k2 and k3","v=",v,"k1=",k1,"k2=",k2,"k3=",k3,"i=",i )
                     print (H<-CBND_3diffsize(k=c(k1,k2,k3),i=i,D=8,sk2=1,sk3=1))
                    }else
                      if((v-(4*k2)-(2*k3)-2)%%(2*k1)==0 & (i=(v-(4*k2)-(2*k3)-2)/(2*k1)) & (i%%1==0) &(i>0)){
                        cat("MCSPBNDS-II is possible with two sets of k2 and one set of k3","v=",v,"k1=",k1,"k2=",k2,"k3=",k3,"i=",i )
                        print(H<-CBND_3diffsize(k=c(k1,k2,k3),i=i,D=8,sk2=2,sk3=1))
                      }else
                        if((v-(2*k2)-(4*k3)-2)%%(2*k1)==0 & (i=(v-(2*k2)-(4*k3)-2)/(2*k1)) & (i%%1==0) &(i>0)){
                          cat("MCSPBNDS-II is possible with one set of k2 and two sets of k3","v=",v,"k1=",k1,"k2=",k2,"k3=",k3,"i=",i )
                          print(H<-CBND_3diffsize(k=c(k1,k2,k3),i=i,D=8,sk2=1,sk3=2))
                        }else
                          if((v-(4*k2)-(4*k3)-2)%%(2*k1)==0 &(i=(v-(4*k2)-(4*k3)-2)/(2*k1)) & (i%%1==0) &(i>0)){
                            cat("MCSPBNDS-II is possible with two sets of k2 and k3","v=",v,"k1=",k1,"k2=",k2,"k3=",k3,"i=",i )
                            print(H<-CBND_3diffsize(k=c(k1,k2,k3),i=i,D=8,sk2=2,sk3=2))
                          }
                            if((v-(2*k2)-(2*k3)-4)%%(2*k1)==0 &(i=(v-(2*k2)-(2*k3)-4)/(2*k1)) & (i%%1==0) &(i>0)){
                              cat("MCPBNDs-II is possible with one sets of k2 and k3","v=",v,"k1=",k1,"k2=",k2,"k3=",k3,"i=",i)
                              print(H<-CBND_3diffsize(k=c(k1,k2,k3),i=i,D=4,sk2=1,sk3=1))
                            }else
                              if((v-(4*k2)-(2*k3)-4)%%(2*k1)==0 &(i=(v-(4*k2)-(2*k3)-4)/(2*k1)) & (i%%1==0) &(i>0)){
                                cat("MCPBNDs-II is possible with two sets of k2 and one set of k3","v=",v,"k1=",k1,"k2=",k2,"k3=",k3,"i=",i )
                                print(H<-CBND_3diffsize(k=c(k1,k2,k3),i=i,D=4,sk2=2,sk3=1))
                              }else
                                if((v-(2*k2)-(4*k3)-4)%%(2*k1)==0 &(i=(v-(2*k2)-(4*k3)-4)/(2*k1)) & (i%%1==0) &(i>0)){
                                  cat("MCPBNDs-II is possible with one set of k2 and two sets of k3","v=",v,"k1=",k1,"k2=",k2,"k3=",k3,"i=",i )
                                  print(H<-CBND_3diffsize(k=c(k1,k2,k3),i=i,D=4,sk2=1,sk3=2))
                                }else
                                  if((v-(4*k2)-(4*k3)-4)%%(2*k1)==0 &(i=(v-(4*k2)-(4*k3)-4)/(2*k1)) & (i%%1==0) &(i>0)){
                                    cat("MCPBNDs-II is possible with two set of k2 and k3","v=",v,"k1=",k1,"k2=",k2,"k3=",k3,"i=",i )
                                    print(H<-CBND_3diffsize(k=c(k1,k2,k3),i=i,D=4,sk2=2,sk3=2))
                                  }
            
            if((v-(2*k2)-(2*k3))%%(2*k1)==0 &(i=(v-(2*k2)-(2*k3))/(2*k1)) & (i%%1==0) &(i>0)& (((v-2)/2)%%4==2 |((v-2)/2)%%4==3)){
              cat("MCWBNDs-I is possible with one set of k2 and k3","v=",v,"k1=",k1,"k2=",k2,"k3=",k3,"i=",i )
              print(H<-CBND_3diffsize(k=c(k1,k2,k3),i=i,D=5,sk2=1,sk3=1))
            }else
              if((v-(4*k2)-(2*k3))%%(2*k1)==0 &(i=(v-(4*k2)-(2*k3))/(2*k1)) & (i%%1==0) &(i>0)& (((v-2)/2)%%4==2 |((v-2)/2)%%4==3)){
                cat("MCWBNDs-I is possible with two sets of k2 and one set of k3","v=",v,"k1=",k1,"k2=",k2,"k3=",k3,"i=",i )
                print(H<-CBND_3diffsize(k=c(k1,k2,k3),i=i,D=5,sk2=2,sk3=1))
              }else
                if((v-(2*k2)-(4*k3))%%(2*k1)==0 &(i=(v-(2*k2)-(4*k3))/(2*k1)) & (i%%1==0) &(i>0)& (((v-2)/2)%%4==2 |((v-2)/2)%%4==3)){
                  cat("MCWBNDs-I is possible with one set of k2 and two sets of k3","v=",v,"k1=",k1,"k2=",k2,"k3=",k3,"i=",i )
                  print(H<-CBND_3diffsize(k=c(k1,k2,k3),i=i,D=5,sk2=1,sk3=2))
                }else
                  if((v-(4*k2)-(4*k3))%%(2*k1)==0 &(i=(v-(2*k2)-(4*k3))/(2*k1)) & (i%%1==0) &(i>0) & (((v-2)/2)%%4==2 |((v-2)/2)%%4==3)){
                    cat("MCWBNDs-I is  possible with two set of k2 and k3","v=",v,"k1=",k1,"k2=",k2,"k3=",k3,"i=",i )
                    print(H<-CBND_3diffsize(k=c(k1,k2,k3),i=i,D=5,sk2=2,sk3=2))
                  }
            
                                    if((v-(2*k2)-(2*k3))%%(2*k1)==0 &(i=(v-(2*k2)-(2*k3))/(2*k1)) & (i%%1==0) &(i>0) & (((v-2)/2)%%4==0 |((v-2)/2)%%4==3)){
                                      cat("MCSPBNDS-I is possible with one set of k2 and k3","v=",v,"k1=",k1,"k2=",k2,"k3=",k3,"i=",i )
                                      print(H<-CBND_3diffsize(k=c(k1,k2,k3),i=i,D=7,sk2=1,sk3=1))
                                    }else
                                      if((v-(4*k2)-(2*k3))%%(2*k1)==0 &(i=(v-(4*k2)-(2*k3))/(2*k1)) & (i%%1==0) &(i>0)& (((v-2)/2)%%4==0 |((v-2)/2)%%4==3)){
                                        cat("MCSPBNDS-I is possible with two sets of k2 and one set of k3","v=",v,"k1=",k1,"k2=",k2,"k3=",k3,"i=",i )
                                        print(H<-CSPBND_3diffsize(k=c(k1,k2,k3),i=i,D=7,sk2=2,sk3=1))
                                      }else
                                        if((v-(2*k2)-(4*k3))%%(2*k1)==0 &(i=(v-(2*k2)-(4*k3))/(2*k1)) & (i%%1==0) &(i>0) & (((v-2)/2)%%4==0 |((v-2)/2)%%4==3)){
                                          cat("MCSPBNDS-I is possible with one set of k2 and two sets of k3","v=",v,"k1=",k1,"k2=",k2,"k3=",k3,"i=",i )
                                          print(H<-CBND_3diffsize(k=c(k1,k2,k3),i=i,D=7,sk2=1,sk3=2))
                                        }else
                                          if((v-(4*k2)-(4*k3))%%(2*k1)==0 &(i=(v-(4*k2)-(4*k3))/(2*k1)) & (i%%1==0) &(i>0) & (((v-2)/2)%%4==0 |((v-2)/2)%%4==3)){
                                            cat("MCSPBNDS-I is possible with two set of k2 and k3","v=",v,"k1=",k1,"k2=",k2,"k3=",k3,"i=",i )
                                            print(H<-CBND_3diffsize(k=c(k1,k2,k3),i=i,D=7,sk2=2,sk3=2))
                                          }
                                              if((v-(2*k2)-(2*k3)+2)%%(2*k1)==0 &(i=(v-(2*k2)-(2*k3)+2)/(2*k1)) & (i%%1==0) &(i>0)& (((v-2)/2)%%4==2 |((v-2)/2)%%4==3)){
                                                              cat("MCSGBNDS-I is possible with one set of k2 and k3","v=",v,"k1=",k1,"k2=",k2,"k3=",k3,"i=",i )
                                                              print(H<-CBND_3diffsize(k=c(k1,k2,k3),i=i,D=9,sk2=1,sk3=1))
                                                            }else
                                                              if((v-(4*k2)-(2*k3)+2)%%(2*k1)==0 &(i=(v-(4*k2)-(2*k3)+2)/(2*k1)) & (i%%1==0) &(i>0) & (((v-2)/2)%%4==2 |((v-2)/2)%%4==3)){
                                                                cat("MCSGBNDS-I is possible with two set of k2 and one set of k3","v=",v,"k1=",k1,"k2=",k2,"k3=",k3,"i=",i )
                                                                print(H<-CBND_3diffsize(k=c(k1,k2,k3),i=i,D=9,sk2=2,sk3=1))
                                                              }else
                                                                if((v-(2*k2)-(4*k3)+2)%%(2*k1)==0 &(i=(v-(2*k2)-(4*k3)+2)/(2*k1)) & (i%%1==0) &(i>0)& (((v-2)/2)%%4==2 |((v-2)/2)%%4==3)){
                                                                  cat("MCSGBNDS-I is possible with one set of k2 and two set of k3","v=",v,"k1=",k1,"k2=",k2,"k3=",k3,"i=",i )
                                                                  print(H<-CBND_3diffsize(k=c(k1,k2,k3),i=i,D=9,sk2=1,sk3=2))
                                                                }else
                                                                  if((v-(4*k2)-(4*k3)+2)%%(2*k1)==0 &(i=(v-(4*k2)-(4*k3)+2)/(2*k1)) & (i%%1==0) &(i>0)& (((v-2)/2)%%4==2 |((v-2)/2)%%4==3)){
                                                                cat("MCSGBNDS-I is possible with two set of k2 and k3","v=",v,"k1=",k1,"k2=",k2,"k3=",k3,"i=",i )
                                                                print(H<-CBND_3diffsize(k=c(k1,k2,k3),i=i,D=9,sk2=2,sk3=2))
                                                             }
                    if((v-(2*k2)-(2*k3)+2)%%(2*k1)==0 &(i=(v-(2*k2)-(2*k3)+2 )/(2*k1)) & (i%%1==0) &(i>0)){
                    cat("MCWBNDs-II is possible with one set of k2 and k3","v=",v,"k1=",k1,"k2=",k2,"k3=",k3,"i=",i )
                    print(H<-CBND_3diffsize(k=c(k1,k2,k3),i=i,D=6,sk2=1,sk3=1))
                    }else
                    if((v-(4*k2)-(2*k3) +2)%%(2*k1)==0 &(i=(v-(4*k2)-(2*k3)+2)/(2*k1)) & (i%%1==0) &(i>0)){
                    cat("MCWBNDs-II is possible with two sets of k2 and one set of k3","v=",v,"k1=",k1,"k2=",k2,"k3=",k3,"i=",i )
                    print(H<-CBND_3diffsize(k=c(k1,k2,k3),i=i,D=6,sk2=2,sk3=1))
                    }else
                    if((v-(2*k2)-(4*k3)+2)%%(2*k1)==0 &(i=(v-(2*k2)-(4*k3)+2)/(2*k1)) & (i%%1==0) &(i>0)){
                    cat("MCWBNDs-II is possible with one set of k2 and two sets of k3","v=",v,"k1=",k1,"k2=",k2,"k3=",k3,"i=",i )
                    print(H<-CBND_3diffsize(k=c(k1,k2,k3),i=i,D=6,sk2=1,sk3=2))
                    }else
                    if((v-(4*k2)-(4*k3) +2)%%(2*k1)==0 &(i=(v-(2*k2)-(4*k3)+2)/(2*k1)) & (i%%1==0) &(i>0)){
                    cat("MCWBNDs-II is  possible with two set of k2 and k3","v=",v,"k1=",k1,"k2=",k2,"k3=",k3,"i=",i )
                    print(H<-CBND_3diffsize(k=c(k1,k2,k3),i=i,D=6,sk2=2,sk3=2))
                    }
                    if((v-(2*k2)-(2*k3)+4)%%(2*k1)==0 &(i=(v-(2*k2)-(2*k3)+4)/(2*k1)) & (i%%1==0) &(i>0) ){
                    cat("MCGBNDs-II is possible with one set of k2 and k3","v=",v,"k1=",k1,"k2=",k2,"k3=",k3,"i=",i )
                    print(H<-CBND_3diffsize(k=c(k1,k2,k3),i=i,D=10,sk2=1,sk3=1))
                    }else
                    if((v-(4*k2)-(2*k3)+4)%%(2*k1)==0 &(i=(v-(4*k2)-(2*k3)+4)/(2*k1)) & (i%%1==0) &(i>0)){
                    cat("MCGBNDs-II is possible with two set of k2 and one set of k3","v=",v,"k1=",k1,"k2=",k2,"k3=",k3,"i=",i )
                    print(H<-CBND_3diffsize(k=c(k1,k2,k3),i=i,D=10,sk2=2,sk3=1))                                                                                                                                  
                    }else
                    if((v-(2*k2)-(4*k3)+4)%%(2*k1)==0 &(i=(v-(2*k2)-(4*k3)+4)/(2*k1)) & (i%%1==0) &(i>0)){
                    cat("MCGBNDs-II is possible with one set of k2 and two sets of k3","v=",v,"k1=",k1,"k2=",k2,"k3=",k3,"i=",i )
                    print(H<-CBND_3diffsize(k=c(k1,k2,k3),i=i,D=10,sk2=1,sk3=2))
                    }else
                     if((v-(4*k2)-(4*k3)+4)%%(2*k1)==0 &(i=(v-(4*k2)-(4*k3)+4)/(2*k1)) & (i%%1==0) &(i>0)){
                    cat("MCGBNDs-II is possible with two set of k2 and k3","v=",v,"k1=",k1,"k2=",k2,"k3=",k3,"i=",i )                                       
                       print(H<-CBND_3diffsize(k=c(k1,k2,k3),i=i,D=10,sk2=2,sk3=2))
                   }
            if((i%%1<0) |(i%%1>0)){ cat("Other designs are not possible","v=",v,"k1=",k1,"k2=",k2,"k3=",k3,"i=",i)
          }}
      }
                                                               (H<-Designs(v=v,k1=k1,k2=k2,k3=k3))
    }
}
#input v,k1>2,k2=0, and k3=0 for equale designs     
#input v,k1>k2,k2>2 and k3=0 for two different designs     
#input v,k1>k2,k2>k3 and k3>2 for three different designs     

(H<-Designs(v=32,k1=6,k2=5,k3=0))



