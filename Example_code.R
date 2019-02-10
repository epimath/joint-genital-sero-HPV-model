# Create data fits from APC results
library("splines")
library("deSolve")
library("Bhat")
library("Epi")
library("ggplot2")

#Data formated for disease model
pop=read.csv("Pop.csv",header=TRUE)[,2:5]
cases1=read.csv("Cases_G-_S-.csv",header=TRUE)[,2:5] #G-,S-
cases2=read.csv("Cases_G+_S-.csv",header=TRUE)[,2:5] #G+,S-
cases3=read.csv("Cases_G-_S+.csv",header=TRUE)[,2:5] #G-,S+
cases4=read.csv("Cases_G+_S+.csv",header=TRUE)[,2:5] #G+,S+
cases13= cases1+cases3 #G-
cases24= cases2+cases4 #G+
cases12= cases1+cases2 #S-
cases34= cases3+cases4 #S+

#Data formated for APC
G_prevdat=as.data.frame(read.csv('APC_G+.csv',header=TRUE))
G_prevdat$C=G_prevdat$P-G_prevdat$A
G_prevdat=G_prevdat[order(G_prevdat$C),]

S_prevdat=as.data.frame(read.csv('APC_S+.csv',header=TRUE))
S_prevdat$C=S_prevdat$P-S_prevdat$A
S_prevdat=S_prevdat[order(S_prevdat$C),]

#############################################################
# APC analyses
#############################################################

#This APC function is adapted from the 'Epi' package and the calcuation 
# of CIs still requres a function from the package. 

APC = function(data,na,nc,ref.c,alpha){
  data <- data[, c("A", "P", "D", "Y")]
  data <- data[complete.cases(data), ]
  A <- data$A
  P <- data$P
  D <- data$D
  Y <- data$Y
  
  c0=ref.c
  
  MA <- ns(A, df = na)
  MA <- cbind(1, MA)
  MC=ns(P-A, df = nc)
  Rc=predict(MC,1980)
  MCr = MC-t(matrix(Rc,(nc), length(data$A)))
  
  A.pt <- unique(A)
  A.pos <- match(A.pt, A)
  C.pt <- unique(P - A)
  C.pos <- match(C.pt, P - A)
  
  
  m.ac= glm(formula = cbind(D, Y - D) ~ -1 + MA + MCr, family = binomial)
  
  
  BCe=function(k,bys) exp(sum(MCr[k,]*bys))
  Ae=function(k,ays) exp(sum(MA[k,]*ays))
  
  a.eff=apply(as.array(1:length(A.pt)),1,Ae,m.ac$coefficients[1:(1+na)])
  Age <- cbind(Age = A.pt, a.eff,ci.exp(m.ac, subset = "MA", 
                                        ctr.mat = MA[A.pos, , drop = FALSE], alpha = alpha))[order(A.pt), ]
  
  c.eff=apply(as.array(C.pos),1,BCe,m.ac$coefficients[(2+na):(1+na+nc)])
  Coh <- cbind(Cph = C.pt, c.eff,ci.exp(m.ac, subset = "MCr", 
                                        ctr.mat = MCr[C.pos, , drop = FALSE], alpha = alpha))[order(C.pt), ]
  
  return(list(Age=Age,Coh=Coh))
}

na=4 #Number of degrees of freedom for age
nc=5 #Number of degrees of freedom for birth cohort
c0=1980 #Reference birth cohort

#Calculate effects
G_prevAC=APC(G_prevdat, na, nc, c0, 0.05)
S_prevAC=APC(S_prevdat, na, nc, c0, 0.05)

#Plot fits to data


DATA=matrix(NA,nrow=378,ncol=2)
DATA[,1]=c(G_prevdat$A,rep(18:59,5))
DATA[1:12,2] = G_prevdat[1:12,3]/G_prevdat[1:12,4]
DATA[13:53,2]= G_prevdat[13:53,3]/G_prevdat[13:53,4]
DATA[54:92,2] = G_prevdat[54:92,3]/G_prevdat[54:92,4]
DATA[93:132,2]= G_prevdat[93:132,3]/G_prevdat[93:132,4]
DATA[133:166,2] = G_prevdat[133:166,3]/G_prevdat[133:166,4]
DATA[167:168,2] = G_prevdat[167:168,3]/G_prevdat[167:168,4]
DATA[169:210,2] = (G_prevAC$Age[,3]*G_prevAC$Coh[G_prevAC$Coh[,1]=="1945",3])/(1+G_prevAC$Age[,3]*G_prevAC$Coh[G_prevAC$Coh[,1]=="1945",3])
DATA[211:252,2] = (G_prevAC$Age[,3]*G_prevAC$Coh[G_prevAC$Coh[,1]=="1955",3])/(1+G_prevAC$Age[,3]*G_prevAC$Coh[G_prevAC$Coh[,1]=="1955",3])
DATA[253:294,2] = (G_prevAC$Age[,3]*G_prevAC$Coh[G_prevAC$Coh[,1]=="1965",3])/(1+G_prevAC$Age[,3]*G_prevAC$Coh[G_prevAC$Coh[,1]=="1965",3])
DATA[295:336,2] = (G_prevAC$Age[,3]*G_prevAC$Coh[G_prevAC$Coh[,1]=="1975",3])/(1+G_prevAC$Age[,3]*G_prevAC$Coh[G_prevAC$Coh[,1]=="1975",3])
DATA[337:378,2] = (G_prevAC$Age[,3]*G_prevAC$Coh[G_prevAC$Coh[,1]=="1985",3])/(1+G_prevAC$Age[,3]*G_prevAC$Coh[G_prevAC$Coh[,1]=="1985",3])
DATA=as.data.frame(DATA)
colnames(DATA)=c("x","y")
DATA$time=factor(c(rep("1944-49",12), rep("1950-59",41), rep("1960-69",39), rep("1970-79",40),
                   rep("1980-89",34),rep("1990-91",2),rep("1944-49",42), rep("1950-59",42), rep("1960-69",42), 
                   rep("1970-79",42), rep("1980-89",42)))
DATA$type=factor(c(rep("Data",168),rep("Model",210)))

colors=c("firebrick2","darkorange","green2","dodgerblue","slateblue","purple2")

q1 = ggplot(data=DATA,aes(x=x,y=y,color=time))+
  geom_point(data=DATA[DATA$type=="Data",],cex=2)+
  scale_color_manual(guide=guide_legend(title="", keyheight=0.45, label.theme=element_text(size=7,angle=0),override.aes = (list(linetype=c(0,0,0,0,0,0)))),values=c(colors,colors[1:5]))+
  geom_line(data=DATA[DATA$type=="Model",],cex=1)+
  theme_bw(base_size = 8) + labs(x='Age', y='Modeled genital prevalence')+
  theme(axis.title.y=element_text(vjust=1.5))+
  theme(legend.position = c(0.8,0.8),legend.key = element_blank(),legend.background=element_blank(),legend.text=element_text(size=7))+
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(color = "black",size=0.4),
        axis.line.y = element_line(color = "black",size=0.4))+
  annotate("text",label="Genital-only age-cohort model",x=42,y=0.4875)+
  scale_x_continuous(limits=c(17,60),expand=c(0,0.001))+
  scale_y_continuous(limits=c(0,0.5),expand=c(0,0.005))#,breaks=c(0,0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16),labels=c('0.0','0.02','0.04','0.06','0.08','0.10','0.12','0.14','0.16'))


print(q1)

DATA=matrix(NA,nrow=378,ncol=2)
DATA[,1]=c(S_prevdat$A,rep(18:59,5))
DATA[1:12,2] = S_prevdat[1:12,3]/S_prevdat[1:12,4]
DATA[13:53,2]= S_prevdat[13:53,3]/S_prevdat[13:53,4]
DATA[54:92,2] = S_prevdat[54:92,3]/S_prevdat[54:92,4]
DATA[93:132,2]= S_prevdat[93:132,3]/S_prevdat[93:132,4]
DATA[133:166,2] = S_prevdat[133:166,3]/S_prevdat[133:166,4]
DATA[167:168,2] = S_prevdat[167:168,3]/S_prevdat[167:168,4]
DATA[169:210,2] = (S_prevAC$Age[,3]*S_prevAC$Coh[S_prevAC$Coh[,1]=="1945",3])/(1+S_prevAC$Age[,3]*S_prevAC$Coh[S_prevAC$Coh[,1]=="1945",3])
DATA[211:252,2] = (S_prevAC$Age[,3]*S_prevAC$Coh[S_prevAC$Coh[,1]=="1955",3])/(1+S_prevAC$Age[,3]*S_prevAC$Coh[S_prevAC$Coh[,1]=="1955",3])
DATA[253:294,2] = (S_prevAC$Age[,3]*S_prevAC$Coh[S_prevAC$Coh[,1]=="1965",3])/(1+S_prevAC$Age[,3]*S_prevAC$Coh[S_prevAC$Coh[,1]=="1965",3])
DATA[295:336,2] = (S_prevAC$Age[,3]*S_prevAC$Coh[S_prevAC$Coh[,1]=="1975",3])/(1+S_prevAC$Age[,3]*S_prevAC$Coh[S_prevAC$Coh[,1]=="1975",3])
DATA[337:378,2] = (S_prevAC$Age[,3]*S_prevAC$Coh[S_prevAC$Coh[,1]=="1985",3])/(1+S_prevAC$Age[,3]*S_prevAC$Coh[S_prevAC$Coh[,1]=="1985",3])
DATA=as.data.frame(DATA)
colnames(DATA)=c("x","y")
DATA$time=factor(c(rep("1944-49",12), rep("1950-59",41), rep("1960-69",39), rep("1970-79",40),
                   rep("1980-89",34),rep("1990-91",2),rep("1944-49",42), rep("1950-59",42), rep("1960-69",42), 
                   rep("1970-79",42), rep("1980-89",42)) )
DATA$type=factor(c(rep("Data",168),rep("Model",210)))
colors=c("firebrick2","darkorange","green2","dodgerblue","slateblue","purple2")



q2 = ggplot(data=DATA,aes(x=x,y=y,color=time))+
  geom_point(data=DATA[DATA$type=="Data",],cex=2)+
  scale_color_manual(guide=guide_legend(title="", keyheight=0.45, label.theme=element_text(size=7,angle=0),override.aes = (list(linetype=c(0,0,0,0,0,0)))),values=c(colors,colors[1:5]))+
  geom_line(data=DATA[DATA$type=="Model",],cex=1)+
  theme_bw(base_size = 8) + labs(x='Age', y='Modeled seroprevalence')+
  theme(axis.title.y=element_text(vjust=1.5))+
  theme(legend.position = c(0.8,0.8),legend.key = element_blank(),legend.background=element_blank(),
        legend.text=element_text(size=7))+
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(color = "black",size=0.4),
        axis.line.y = element_line(color = "black",size=0.4))+
  annotate("text",label="Sero-only age-cohort model",x=43,y=0.97)+
  scale_x_continuous(limits=c(17,60),expand=c(0,0.001))+
  scale_y_continuous(limits=c(0,1),expand=c(0,0.005))#,breaks=c(0,0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16),labels=c('0.0','0.02','0.04','0.06','0.08','0.10','0.12','0.14','0.16'))


print(q2)

#############################################################
# Partner acquisition rate (adapted from Ryser et al., 2017)
#############################################################
prev=c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
       0, 0, 0, 0, 0, 0, 0.12375, 0.2475, 0.37125, 0.495, 0.403571429, 
       0.312142857, 0.220714286, 0.1637, 0.1411, 0.1185, 0.0959, 0.0733, 
       0.0601, 0.0563, 0.0525, 0.0487, 0.0449, 0.0434, 0.0442, 0.045, 
       0.0458, 0.0466, 0.0488, 0.0524, 0.056, 0.0596, 0.0632, 0.061, 
       0.053, 0.045, 0.037, 0.029, 0.024, 0.0226, 0.021, 0.0194, 0.0178, 
       0.0186, 0.0218, 0.025, 0.0282, 0.0314, 0.0346, 0.0378)

DATA=matrix(NA,nrow=60,ncol=4)
DATA[,1]=0:59
DATA[,2] = prev
DATA=as.data.frame(DATA)
colnames(DATA)=c("x","y","ymin","ymax")

p = ggplot(data=DATA,aes(x=x,y=y))+
  geom_line(cex=0.75)+
  scale_color_manual(guide=guide_legend(title="",label.theme=element_text(size=7,angle=0)),values=c("gray50"))+
  theme_bw(base_size = 8) + labs(x='Age', y='Partner acquisition rate')+
  theme(axis.title.y=element_text(vjust=1.5))+
  theme(legend.position = c(0.75,0.90),legend.key = element_blank(),legend.background=element_blank(),legend.text=element_text(size=7))+
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(color = "black",size=0.4),
        axis.line.y = element_line(color = "black",size=0.4))+
  scale_x_continuous()+
  scale_y_continuous()

print(p)


##############################################################
## Disease model
##############################################################

#Linear interpolation function will be used
#to get continuous partner acquisition rate as a function of age
lin_interp=function(u,P){  (P[(ceiling(u)+1)]-P[(floor(u)+1)])*(u-floor(u))+P[(floor(u)+1)]}

#Differential equation
SI_sero=function(t,x,par)
{
  lambda0=par[1]   #Force of infection at reference value
  rho=par[2]       #Reduction in probability of infection when seropositive
  gamma_neg=par[3] #Recovery rate when seronegative
  gamma_pos=par[4] #Recovery rate when seropositive
  sigma=par[5]     #Serocoersion rate
  omega=par[6]     #Antibody decay rate
  nu=par[7]        #Rate of entering latency
  mu=par[8]        #Rate of reactivation
  lambdat=lambda0*lin_interp(t,prev) #Force of infection at age t
  S_neg = x[1]
  I_neg = x[2]
  L_neg = x[3]
  S_pos = x[4]
  I_pos = x[5]
  L_pos = x[6]
  dSnegdt= gamma_neg*I_neg + omega*S_pos - lambdat*S_neg
  dInegdt= lambdat*S_neg + omega*I_pos + mu*L_neg - (gamma_neg + sigma + nu)*I_neg 
  dLnegdt= nu*I_neg + omega*L_pos - mu*L_neg
  dSposdt= gamma_pos*I_pos -(rho*lambdat + omega)*S_pos
  dIposdt= rho*lambdat*S_pos + sigma*I_neg + mu*L_pos - (gamma_pos + omega + nu)*I_pos
  dLposdt= nu*I_pos - (mu + omega)*L_pos
  
  
  dx=list(c(dSnegdt,dInegdt,dLnegdt,dSposdt, dIposdt, dLposdt))
  return(dx)
}

############################################################
### Fit model for all cohorts using real data and splines
############################################################
#Set up splines for prevalence
my_ns=function(x,df,intercept=FALSE){
  n=length(x)
  t=seq(from=x[1],to=x[n],by=((n-1)/df)) #knots
  k=length(t)  #number of knots
  intercept=FALSE  #include an intercept?
  X=matrix(1,n,k-1L+intercept)
  scale= (t[k] - t[1])^(2/3)#scaling factor to avoid ill-conditioning. See rcspline.eval
  X[,1L+intercept]=x/scale
  for (j in 1:(k-2L)){
    X[,(j+1+intercept)]=pmax((x-t[j])/scale,0)^3- pmax((x-t[k-1])/scale,0)^3*(t[k]-t[j])/(t[k]-t[k-1]) + pmax((x-t[k])/scale,0)^3*(t[k-1]-t[j])/(t[k]-t[k-1])  
  }
  rownames(X)=x
  
  return(X)
}

nc = 5 #Number of degrees of freedom. Approximately one per ten years is a good starting point.
MC=my_ns(1944:1991,nc)
Rcs=MC["1980",] #reference cohort

#Function for returning age-specific effects from the spline parameters
BCe=function(k,bcps) exp(sum((MC[k,]-Rcs)*bcps))
#Bounds on spline parameters for optimization step
bclab=paste("BCspline", 1:nc, sep = "")
bcl=rep(-10,nc)
bcu=rep(10,nc)

#Multinomial negative log-likelihood
NLL_joint=function(par){
  lambda0=par[1]
  rho=1
  gamma_neg=par[2]
  gamma_pos=par[2]
  sigma=par[3]
  omega=par[4]
  nu=par[5]
  mu=par[6]
  bcps=par[7:(6+nc)]
  bceff=apply(as.array(1:length(1944:1991)),1,BCe,bcps)
  
  #set up matrices for modeled prevalences
  P1=matrix(NA,length(0:59),length(1944:1991))
  P2=matrix(NA,length(0:59),length(1944:1991))
  P3=matrix(NA,length(0:59),length(1944:1991))
  P4=matrix(NA,length(0:59),length(1944:1991))
  #For each cohort, calculate hazard based on prevalence
  for (i in 1:length(1944:1991)){
    lambda=lambda0*bceff[i]
    out=ode(y=c(Y1=1,Y2=0,Y3=0,Y4=0,Y5=0,Y6=0),seq(0,59,1),SI_sero, c(lambda, rho, gamma_neg, gamma_pos, sigma,omega, nu,mu), method='lsode',atol=1e-15) 
    P1[,i]=out[,2]+out[,4]
    P2[,i]=out[,3]
    P3[,i]=out[,5]+out[,7]
    P4[,i]=out[,6]
  }
  
  
  H1 = matrix(0,42,7)
  #in H1: age = i+17, period = j + 2002, k=cohort-1943 = (j+2002)-(i+17)-1943=j-i+42
  for (i in 1:42){
    for (j in 1:7){
      if ( ((j-i)>=-41)){
        H1[i,j]=P1[i+18,j-i+42]
      }}}
  H2 = matrix(0,42,7)
  for (i in 1:42){
    for (j in 1:7){
      if ( ((j-i)>=-41)){
        H2[i,j]=P2[i+18,j-i+42]
      }}}
  H3 = matrix(0,42,7)
  for (i in 1:42){
    for (j in 1:7){
      if ( ((j-i)>=-41)){
        H3[i,j]=P3[i+18,j-i+42]
      }}}
  H4 = matrix(0,42,7)
  for (i in 1:42){
    for (j in 1:7){
      if ( ((j-i)>=-41)){
        H4[i,j]=P4[i+18,j-i+42]
      }}}
  
  fNLL=-sum(cases1*log(H1[,c(1,3,5,7)])+cases2*log(H2[,c(1,3,5,7)])+cases3*log(H3[,c(1,3,5,7)])+cases4*log(H4[,c(1,3,5,7)]))
  
  return(fNLL)
}


#Optimization step
x = list(label=c('lambda0','gamma','simga','omega','nu','mu',bclab),
          est=c(0.51, 0.41, 0.74, 0.047, 1.06, 0.53, 0.25, -0.014, -0.085, 0.12, 0.85),
          low=c(rep(1e-3,6),bcl),upp=c(rep(1E1,6),bcu))
NLL_joint(x$est)
#mle= dfp(x,NLL_joint) #commented out for ease of code running

#Output of optimization step - your exact solution may vary depending on initial estimate and your computer
#est_j=mle$est
est_j=c(0.508054164365755, 0.405722640026515, 0.742423072097424, 0.0479606229421793, 
        1.05713960979968, 0.526112711932872, 0.249642353611332, -0.0141234424837452, 
        -0.0852690987459609, 0.116225604732557, 0.848525490546582)
NLL_joint(est_j)


#Solve for Hessian to get confidence bounds
#out_joint=optim(est_j,NLL_joint,hessian=TRUE)
out_joint = structure(c(5369.21658499523, -3803.12448078257, 2014.35262488303, 
            -30602.4553487987, -298.50968962819, 591.501680560214, -2065.95506654139, 
            -23018.4500427413, -9635.01205603734, -2353.7273450529, 1.33348430608748, 
            -3803.12448078257, 3857.38060413132, -1191.12413881339, 22189.5992226564, 
            236.224254194894, -437.645922488628, 1966.42439937023, 23725.9481003775, 
            10338.7401422879, 2855.72342500018, 205.141594278757, 2014.35262488303, 
            -1191.12413881339, 1206.59810818324, -14167.0133812113, -159.188979750979, 
            254.240721460519, -693.492747018354, -7858.95135675219, -3261.65966737335, 
            -806.361650347753, 32.4323688118966, -30602.4553487987, 22189.5992226564, 
            -14167.0133812113, 295235.123122893, -34.7665765048077, -726.384994209184, 
            19809.6779880643, 235508.671223101, 103062.221569303, 29436.6298400064, 
            2652.770603504, -298.50968962819, 236.224254194894, -159.188979750979, 
            -34.7665765048077, 158.435752155128, -170.38735825281, -133.368926867661, 
            -1764.78452169704, -899.056720527369, -364.211306191464, -102.106490430742, 
            591.501680560214, -437.645922488628, 254.240721460519, -726.384994209184, 
            -170.38735825281, 321.615371376538, 120.552952921571, 1800.41916974005, 
            901.034627077024, 398.214588699375, 116.9869133264, -2065.95506654139, 
            1966.42439937023, -693.492747018354, 19809.6779880643, -133.368926867661, 
            120.552952921571, 1843.23899111405, 21956.2342200561, 9724.76531296707, 
            2903.49468320983, 345.664040992233, -23018.4500427413, 23725.9481003775, 
            -7858.95135675219, 235508.671223101, -1764.78452169704, 1800.41916974005, 
            21956.2342200561, 285759.618875773, 130812.701189143, 41597.1399978616, 
            5883.99818786911, -9635.01205603734, 10338.7401422879, -3261.65966737335, 
            103062.221569303, -899.056720527369, 901.034627077024, 9724.76531296707, 
            130812.701189143, 60769.7641376035, 19837.5444128942, 2993.75108852473, 
            -2353.7273450529, 2855.72342500018, -806.361650347753, 29436.6298400064, 
            -364.211306191464, 398.214588699375, 2903.49468320983, 41597.1399978616, 
            19837.5444128942, 6824.35106000412, 1179.60871943978, 1.33348430608748, 
            205.141594278757, 32.4323688118966, 2652.770603504, -102.106490430742, 
            116.9869133264, 345.664040992233, 5883.99818786911, 2993.75108852473, 
            1179.60871943978, 250.330381959429), .Dim = c(11L, 11L))

cov_joint=solve(out_joint)

#Confidence bounds for splines
low_j=rep(NA,length(1944:1991))
high_j=rep(NA,length(1944:1991))
tMC=t(MC)
for (i in 1:length(1944:1991)){
  low_j[i] = exp(sum((MC[i,]-Rcs)*est_j[7:11]) -qnorm(0.975) *sqrt((MC[i,]-Rcs)%*% cov_joint[7:11,7:11] %*% (tMC[,i]-Rcs)))
  high_j[i] = exp(sum((MC[i,]-Rcs)*est_j[7:11]) +qnorm(0.975) *sqrt((MC[i,]-Rcs)%*% cov_joint[7:11,7:11] %*% (tMC[,i]-Rcs)))
}

#Confidence bounds for model parameters
est_j[1:6]- sqrt(diag(cov_joint)[1:6])*qnorm(0.975)
est_j[1:6]+ sqrt(diag(cov_joint)[1:6])*qnorm(0.975)



# Plot fits to data

lambda0=est_j[1]
rho=1
gamma_neg=est_j[2]
gamma_pos=est_j[2]
sigma=est_j[3]
omega=est_j[4]
nu=est_j[5]
mu=est_j[6]
bcps=est_j[7:(6+nc)]
bceff_joint=apply(as.array(1:length(1944:1991)),1,BCe,bcps)

#set up matrices for modeled prevalences
P1=matrix(NA,length(0:59),length(1944:1991))
P2=matrix(NA,length(0:59),length(1944:1991))
P3=matrix(NA,length(0:59),length(1944:1991))
P4=matrix(NA,length(0:59),length(1944:1991))

#For each cohort, calculate hazard based on prevalence
for (i in 1:length(1944:1991)){
  lambda=lambda0*bceff_joint[i]
  out=ode(y=c(Y1=1,Y2=0,Y3=0,Y4=0,Y5=0,Y6=0),seq(0,59,1),SI_sero, c(lambda, rho, gamma_neg, gamma_pos, sigma,omega, nu,mu), method='lsode',atol=1e-15) 
  P1[,i]=out[,2]+out[,4]
  P2[,i]=out[,3]
  P3[,i]=out[,5]+out[,7]
  P4[,i]=out[,6]
}

DATA=matrix(NA,nrow=378,ncol=2)
DATA[,1]=c(G_prevdat$A,rep(18:59,5))
DATA[1:12,2] = G_prevdat[1:12,3]/G_prevdat[1:12,4]
DATA[13:53,2]= G_prevdat[13:53,3]/G_prevdat[13:53,4]
DATA[54:92,2] = G_prevdat[54:92,3]/G_prevdat[54:92,4]
DATA[93:132,2]= G_prevdat[93:132,3]/G_prevdat[93:132,4]
DATA[133:166,2] = G_prevdat[133:166,3]/G_prevdat[133:166,4]
DATA[167:168,2] = G_prevdat[167:168,3]/G_prevdat[167:168,4]
DATA[169:210,2] = P2[19:60,2]+P4[19:60,2]
DATA[211:252,2] = P2[19:60,12]+P4[19:60,12]
DATA[253:294,2] = P2[19:60,22]+P4[19:60,22]
DATA[295:336,2] = P2[19:60,32]+P4[19:60,32]
DATA[337:378,2] = P2[19:60,42]+P4[19:60,42]
DATA=as.data.frame(DATA)
colnames(DATA)=c("x","y")
DATA$time=factor(c(rep("1944-49",12), rep("1950-59",41), rep("1960-69",39), rep("1970-79",40),
                   rep("1980-89",34),rep("1990-91",2),rep("1944-49",42), rep("1950-59",42), rep("1960-69",42), 
                   rep("1970-79",42), rep("1980-89",42)) )
DATA$type=factor(c(rep("Data",168),rep("Model",210)))

colors=c("firebrick2","darkorange","green2","dodgerblue","slateblue","purple2")

p1 = ggplot(data=DATA,aes(x=x,y=y,color=time))+
  geom_point(data=DATA[DATA$type=="Data",],cex=2)+
  scale_color_manual(guide=guide_legend(title="", keyheight=0.45, label.theme=element_text(size=7,angle=0),override.aes = (list(linetype=c(0,0,0,0,0,0)))),values=c(colors,colors[1:5]))+
  geom_line(data=DATA[DATA$type=="Model",],cex=1)+
  theme_bw(base_size = 8) + labs(x='Age', y='Modeled genital prevalence')+
  theme(axis.title.y=element_text(vjust=1.5))+
  theme(legend.position = c(0.8,0.8),legend.key = element_blank(),legend.background=element_blank(),legend.text=element_text(size=7))+
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(color = "black",size=0.4),
        axis.line.y = element_line(color = "black",size=0.4))+
  annotate("text",label="Joint genital-sero disease model",x=41,y=0.4875)+
  scale_x_continuous(limits=c(17,60),expand=c(0,0.001))+
  scale_y_continuous(limits=c(0,0.5),expand=c(0,0.005))

print(p1)

DATA=matrix(NA,nrow=378,ncol=2)
DATA[,1]=c(S_prevdat$A,rep(18:59,5))
DATA[1:12,2] = S_prevdat[1:12,3]/S_prevdat[1:12,4]
DATA[13:53,2]= S_prevdat[13:53,3]/S_prevdat[13:53,4]
DATA[54:92,2] = S_prevdat[54:92,3]/S_prevdat[54:92,4]
DATA[93:132,2]= S_prevdat[93:132,3]/S_prevdat[93:132,4]
DATA[133:166,2] = S_prevdat[133:166,3]/S_prevdat[133:166,4]
DATA[167:168,2] = S_prevdat[167:168,3]/S_prevdat[167:168,4]
DATA[169:210,2] = P4[19:60,2]+P3[19:60,2]
DATA[211:252,2] = P4[19:60,12]+P3[19:60,12]
DATA[253:294,2] = P4[19:60,22]+P3[19:60,22]
DATA[295:336,2] = P4[19:60,32]+P3[19:60,32]
DATA[337:378,2] = P4[19:60,42]+P3[19:60,42]
DATA=as.data.frame(DATA)
colnames(DATA)=c("x","y")
DATA$time=factor(c(rep("1944-49",12), rep("1950-59",41), rep("1960-69",39), rep("1970-79",40),
                   rep("1980-89",34),rep("1990-91",2),rep("1944-49",42), rep("1950-59",42), rep("1960-69",42), 
                   rep("1970-79",42), rep("1980-89",42)) ) 
DATA$type=factor(c(rep("Data",168),rep("Model",210)))
colors=c("firebrick2","darkorange","green2","dodgerblue","slateblue","purple2")

p2 = ggplot(data=DATA,aes(x=x,y=y,color=time))+
  geom_point(data=DATA[DATA$type=="Data",],cex=2)+
  scale_color_manual(guide=guide_legend(title="", keyheight=0.45, label.theme=element_text(size=7,angle=0),override.aes = (list(linetype=c(0,0,0,0,0,0)))),values=c(colors,colors[1:5]))+
  geom_line(data=DATA[DATA$type=="Model",],cex=1)+
  theme_bw(base_size = 8) + labs(x='Age', y='Modeled seroprevalence')+
  theme(axis.title.y=element_text(vjust=1.5))+
  theme(legend.position = c(0.8,0.8),legend.key = element_blank(),legend.background=element_blank(),legend.text=element_text(size=7))+
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(color = "black",size=0.4),
        axis.line.y = element_line(color = "black",size=0.4))+
  annotate("text",label="Joint genital-sero disease model",x=41,y=2*0.4875)+
  scale_x_continuous(limits=c(17,60),expand=c(0,0.001))+
  scale_y_continuous(limits=c(0,1),expand=c(0,0.005))#,breaks=c(0,0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16),labels=c('0.0','0.02','0.04','0.06','0.08','0.10','0.12','0.14','0.16'))

print(p2)
