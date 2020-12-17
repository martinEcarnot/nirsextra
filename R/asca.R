asca <-function(...)  {

library(MASS)

  X=dat$x
  desmatrix=fact

  #########################
  gendmat <-function(designmat) {

    ns<-nrow(designmat)
    nfact<-ncol(designmat)
    dmain<- list() #data.frame(nrow=1,ncol=nfact)  # dmain<-cell(1,nfact)

    for (i in 1:nfact) {
      lev<-unique(designmat[,i])
      nl<-length(lev)
      dmat<-matrix(0,ns,nl-1)
      for (j in 1 : nl-1) {
        dmat[designmat[,i]==lev[j],j] <- 1
      } #
      dmat[designmat[,i]==lev[nl],]<--1
      dmain[[i]]<-dmat
    } #
    return(dmain)
  }

####################
  createdesign <-function(designmat) {

    ns<-nrow(designmat)
    nfact<-ncol(designmat)
    nfact<-length(designmat)
    ns<-nrow(designmat[[1]])

    indmat<-expand.grid(replicate(nfact, 1:2, simplify=FALSE))
    nmat<-nrow(indmat)
    dmatrices<- vector("list", nmat)
    desterms<-vector("list", nmat)
    deslabels<-vector("list", nmat)
    desorder<-matrix(0,nrow=nmat,ncol=1)

    for ( i in 1 : nmat) {
      dm<-matrix(1,nrow=ns,ncol=1)
      for ( j in 1 : nfact) {
        if (indmat[i,j]==1) {
          effmat<-designmat[[j]]
        }else{
          effmat<-matrix(1,nrow=ns,ncol=1)
        }
        dm<-kronecker(dm,effmat)
        dm<-dm[seq(1,nrow(dm),ns+1),]
      }
      dmatrices[[i]]<-dm
      desterms[[i]]<-which(indmat[i,]==1)
      deslabels[[i]]<-paste0(LETTERS[which(indmat[i,]==1)],collapse ="")
      desorder[i]<-length(desterms[[i]])
    }

    #Sorting according to increasing order of interactions
    newindex=order(nchar(unlist(deslabels)), unlist(deslabels))
    deslabels = unlist(deslabels)[newindex]

    desterms<-desterms[newindex]
    dmatrices<-dmatrices[newindex]
    desorder<-desorder[newindex]
    # deslabels<-t(cellstr(deslabels))
    deslabels[1]<-'Mean'

    return(list(dmatrices=dmatrices,desterms=desterms, deslabels=deslabels, desorder=desorder))
  }

  ############################"

  args = list(...)
  # nargin <- length(as.list(match.call())) -1

if (nargin==1 && strcmp(args[[1]], 'options')) {
    opt=list()
    opt$preproc<-'none'
    opt$reducedmodel<-'standard'
    opt$permtest<-'on'
    opt$permfacts<-'all'
    opt$nperm<-10000
    opt$bootstrap<-'off'
    opt$bootmatrix<-'original'
    opt$bootsave<-'confint'
    opt$nboot<-1000
    opt$plots<-'on'
    opt$confl<-0.95
    model<-opt
    return(model)
} else if (nargin==2) {
    # X<-varargin[[1]]
    # desmatrix<-varargin{2}
    opt$preproc<-'none'
    opt$reducedmodel<-'standard'
    opt$permtest<-'on'
    opt$permfacts<-'all'
    opt$nperm<-10000
    opt$bootstrap<-'off'
    opt$bootmatrix<-'original'
    opt$bootsave<-'confint'
    opt$nboot<-1000
    opt$plots<-'on'
    opt$confl<-0.95
} else if (nargin==3) {
    # X<-varargin[[1]]
    # desmatrix<-varargin{2}
    # opt<-varargin{3}
 } #


#preprocessing
ntot<-nrow(X)

switch(opt$preproc,
    none={Xp<-X},
    auto={Xp<-X./repmat(std(X), ntot, 1)},
    pareto={Xp<-X./repmat(sqrt(std(X)),ntot, 1)}
)

Xdata=list(PreprData=I(Xp),TotalSSQ=sum(sum(Xp^2)))
model=list(Xdata=Xdata,Design=desmatrix)

dmain<-gendmat(desmatrix)
ldmatrices = createdesign(dmain)
dmatrices=ldmatrices$dmatrices
desterms=ldmatrices$desterms
deslabels=ldmatrices$deslabels
desorder=ldmatrices$desorder

Xd<-Xp

for ( i in 1 : length(dmatrices)) {
    Xeff<-dmatrices[[i]] %*% ginv(dmatrices[[i]]) %*%Xd
    ssqEff<-sum(sum(Xeff^2))
    Xd<-Xd-Xeff
    l<-deslabels[[i]]
    eval(parse(text=paste0('model$X', l,'$EffectMatrix<-Xeff')))
    eval(parse(text=paste0('model$X', l,'$EffectSSQ<-ssqEff')))

    if (i==1) {
        ssqtot<-sum(sum(Xd^2))
        model$Xdata$CenteredData<-Xd
        model$Xdata$CenteredSSQ<-ssqtot
    } else {
        expVar<-100*(ssqEff/ssqtot)
        eval(parse(text=paste0('model$X', l,'$EffectExplVar<-expVar')))
     } #

    eval(parse(text=paste0('model$X', l,'$DesignMatrix<-dmatrices[[i]]')))
    eval(parse(text=paste0('model$X', l,'$DesignTerms<-desterms[[i]]')))
    eval(parse(text=paste0('model$X', l,'$EffectLabel<-deslabels[[i]]')))
    eval(parse(text=paste0('model$X', l,'$TermOrder<-desorder[i]')))

 } #

model$XRes$EffectMatrix<-Xd
model$XRes$EffectSSQ<-sum(sum(Xd^2))
model$XRes$EffectExplVar<-100*(model$XRes$EffectSSQ/ssqtot)
model$XRes$EffectLabel<-'Res'
model$XRes$TermOrder<-max(desorder)+1

for ( i in 2 : length(dmatrices)) {
    l<-deslabels[[i]]

    if (is.character(opt$reducedmodel) && opt$reducedmodel=='standard') {
        if (i==2) {
          remfact<-desterms[(i+1):length(desterms)]
        } else if (i==length(dmatrices)) {
          remfact<-desterms[2:(length(desterms)-1)]
        } else {
          remfact<-desterms[c(2:(i-1),(i+1):length(desterms))]
        }
    } else {
        remfact<-opt$reducedmodel[[i]]
     } #
    Xx<-model$Xdata$CenteredData
    for ( j in 1 : length(remfact)) {
        m<-LETTERS[remfact[[j]]]
        eval(parse(text=paste0('Xx<-Xx-model$X',m,'$EffectMatrix')))
     } #
    eval(parse(text=paste0('model$X', l,'$ReducedMatrix<-Xx')))
 } #

model$TermLabels<-deslabels
model$Options<-opt


#Permutation tests, if required
model<-ascaptest(model)

#SCA Models
model<-ascasca(model)

#bootstrap, if required
model<-ascaboot(model)


}







#
#
#
#
# #############################################
ascaptest <-function(ascamodel) {


pmodel<-ascamodel
dlab<-ascamodel$TermLabels

if  (ascamodel$Options$permtest=='on') {      #  strcmp(ascamodel.Options.permtest, 'on') {

    signfacts<-list(length(dlab)-1)
    sc<-0
}else{
    signfacts<-NULL
 } #


for ( i in 2 : length(dlab)) {
    l<-dlab[[i]]
    if (ascamodel$Options$permtest == 'on') {
        if (is.character(ascamodel$Options$permfacts) && ascamodel$Options$permfacts== 'all') {

            eval(parse(text=paste0('Xr<-ascamodel$X', l, '$ReducedMatrix ')))
            eval(parse(text=paste0('Dr<-ascamodel$X', l, '$DesignMatrix ')))
            ssqp<-ptest(Xr,Dr, ascamodel$Options$nperm)
            eval(parse(text=paste0('seff<-ascamodel$X', l, '$EffectSSQ ')))
            p<-length(which(ssqp>=seff))/ascamodel$Options$nperm
            if (p<<-0.05) {
                sc<-sc+1
                signfacts[[sc]]<-l
             } #


            eval(parse(text=paste0('pmodel$X',l,'$EffectSignif$NullDistr<-ssqp')))
            eval(parse(text=paste0('pmodel$X',l,'$EffectSignif$p<-p')))
        }else{
            if (is.element(l, ascamodel$Options$permfacts)) {  # Or  ascamodel$Options$permfacts == 'all'
                eval(parse(text=paste0('Xr<-ascamodel$X', l, '$ReducedMatrix ')))
                eval(parse(text=paste0('Dr<-ascamodel$X', l, '$DesignMatrix ')))
                ssqp<-ptest(Xr,Dr, ascamodel$Options$nperm)
                eval(parse(text=paste0('seff<-ascamodel$X', l, '$EffectSSQ ')))
                p<-length(find(ssqp>=seff))/ascamodel$Options$nperm
                if (p<<-0.05) {
                    sc<-sc+1
                    signfacts[[sc]]<-l
                 } #

                eval(parse(text=paste0('pmodel$X',l,'$EffectSignif$NullDistr<-ssqp')))
                eval(parse(text=paste0('pmodel$X',l,'$EffectSignif$p<-p')))
            }else{
                eval(parse(text=paste0('pmodel$X',l,'$EffectSignif$NullDistr<-NULL]')))
                eval(parse(text=paste0('pmodel$X',l,'$EffectSignif$p<-NULL')))


             } #
         } #
    }else{
        eval(parse(text=paste0('pmodel$X',l,'$EffectSignif$NullDistr<-[]')))
        eval(parse(text=paste0('pmodel$X',l,'$EffectSignif$p<-[]')))
     } #
 } #

if (ascamodel$Options$permtest == 'on') {
    signfacts<-signfacts[[1:sc]]
 } #

pmodel$SignificantTerms<-signfacts
return(pmodel)
}

# ########################################################
ptest <-function(X,D, nperm) {

ns<-nrow(X)  #Number of samples
ssqp<-matrix(0,nrow=nperm,ncol=1)   #Initialization of the permuted SSQ vector
if (is.null(dim(D))) {dim(D)=c(length(D),1)}

for ( i in 1 : nperm) {
    hh<-sample(ns)
    Xpp<-D[hh,] %*% ginv(D[hh,]) %*% X
    ssqp[i]<-sum(Xpp^2)
 } #

return(ssqp)
}
#
# ##########################################################
ascasca <-function(ascamodel) {

smodel<-ascamodel
dlab<-ascamodel$TermLabels

for ( i in 2 : length(dlab)) {
    l<-dlab[[i]]
    eval(parse(text=paste0('Xr<-ascamodel$X', l, '$EffectMatrix ')))
    eval(parse(text=paste0('ssqr<-ascamodel$X', l, '$EffectSSQ ')))

    R<-qr(Xr)$rank
    rsvd<-svd(Xr,R)
    t<-rsvd$u %*% rsvd$s
    taug<-(Xr+ascamodel$XRes$EffectMatrix) %*% rsvd$v
    varex<-100*(diag(s)^2)/ssqr

    eval(parse(text=paste0('smodel$X', l, '$SCA$Model$Scores<-t ')))
    eval(parse(text=paste0('smodel$X', l, '$SCA$Model$ScoreswithRes<-taug ')))
    eval(parse(text=paste0('smodel$X', l, '$SCA$Model$Loadings<-P ')))
    eval(parse(text=paste0('smodel$X', l, '$SCA$Model$ExplVar<-varex ')))
#       eval(parse(text=paste0('smodel$X', l, '$SCA$Model$s<-s ']) # Only if we want to
#       trace error : Q residuals etc...
 } #

return(smodel)

}

# ##########################################################
ascaboot <-function(ascamodel) {

bmodel<-ascamodel
dlab<-ascamodel$TermLabels
Xd<-ascamodel$Xdata$CenteredData

for ( i in 2 : length(dlab)) {
    l<-dlab[[i]]

    if (ascamodel$Options$bootstrap == 'all') {
        if (ascamodel$Options$bootmatrix == 'original') {

            Xd<-ascamodel$Xdata$CenteredData
        }else if (ascamodel$Options$bootmatrix == 'reduced') {
            eval(parse(text=paste0('Xd<-ascamodel$X', l, '$ReducedMatrix ')))
         } #

        eval(parse(text=paste0('Dd<-ascamodel$X', l, '$DesignMatrix ')))
        eval(parse(text=paste0('Pd<-ascamodel$X', l, '$SCA$Model$Loadings ')))
        resBootload<-bootload(Xd,Dd, Pd, ascamodel$Options$confl, ascamodel$Options$nboot)
        Pb=resBootload$Pb
        Pbcrit=resBootload$Pbcrit
        svars=resBootload$svars

    } else if (ascamodel$Options$bootstrap == 'signif') {
        if (is.element(l, ascamodel$SignificantTerms)) {
            if (ascamodel$Options$bootmatrix == 'original') {

                Xd<-ascamodel$Xdata$CenteredData
            } else if (ascamodel$Options$bootmatrix == 'reduced') {

                eval(parse(text=paste0('Xd<-ascamodel$X', l, '$ReducedMatrix ')))
             }
            eval(parse(text=paste0('Dd<-ascamodel$X', l, '$DesignMatrix ')))
            eval(parse(text=paste0('Pd<-ascamodel$X', l, '$SCA$Model$Loadings ')))
            resBootload<-bootload(Xd,Dd, Pd, ascamodel$Options$confl, ascamodel$Options$nboot)
            Pb=resBootload$Pb
            Pbcrit=resBootload$Pbcrit
            svars=resBootload$svars
        } else {
            Pb<-NULL
            Pbcrit<-NULL
            svars<-NULL
         } #


    } else if (ascamodel$Options$bootstrap == 'off') {
        Pb<-NULL
        Pbcrit<-NULL
        svars<-NULL
     } #



    switch(opt$preproc,
           none={Xp<-X},
           auto={Xp<-X./repmat(std(X), ntot, 1)},
           pareto={Xp<-X./repmat(sqrt(std(X)),ntot, 1)}
    )
    switch (ascamodel$Options$bootsave,
        all={
            eval(parse(text=paste0('bmodel$X', l, '$SCA$Bootstrap$Loadings<-Pb ')))
            eval(parse(text=paste0('bmodel$X', l, '$SCA$Bootstrap$ConfIntervals<-Pbcrit ')))
            eval(parse(text=paste0('bmodel$X', l, '$SCA$Bootstrap$SignificantVariables<-svars ')))
},
        confint={
            eval(parse(text=paste0('bmodel$X', l, '$SCA$Bootstrap$Loadings<-[] ')))
            eval(parse(text=paste0('bmodel$X', l, '$SCA$Bootstrap$ConfIntervals<-Pbcrit ')))
            eval(parse(text=paste0('bmodel$X', l, '$SCA$Bootstrap$SignificantVariables<-svars ')))
},
        signvars={
            eval(parse(text=paste0('bmodel$X', l, '$SCA$Bootstrap$Loadings<-[] ')))
            eval(parse(text=paste0('bmodel$X', l, '$SCA$Bootstrap$ConfIntervals<-[] ')))
            eval(parse(text=paste0('bmodel$X', l, '$SCA$Bootstrap$SignificantVariables<-svars ')))
}
     ) #

 } #
return(bmodel)
}


# #########################################################
bootload <-function(Xd,Dd, Pd, confl, nboot) {
#Bootstrap

sl<-(confl+1)/2
ll<-1-sl
svars<-list(ncol(Pd))

Pb<- array(0, c(nboot, nrow(Pd), ncol(Pd))) # Pb=zeros([nboot size(Pd)]);

bootp<-matrix(0,nrow(Xd))
lev<-unique(Dd, 'rows')

for ( i in 1 : nboot) {
    for ( j in 1 : nrow(lev,1)) {
        xx<-which(is.element(Dd, lev[j,], 'rows')==1)
        yy<-ceiling(length(xx)*sample(length(xx)))  # rand(1,length(xx)))
        bootp(xx)<-xx[yy]
     } #
    Xboot<-Xd[bootp,]
    Xdp<-Dd %*% ginv(Dd) %*% Xboot
    rsvd<-svd(Xdp,ncol(Pd))  # [~,~,vb]<-svds(Xdp,size(Pd,2))
    resOrth.proc<-orth.proc(Pd,rsvd$v)
    Pb[i,,]=resOrth.proc$yrot
 } #
Pb<-sort(Pb)
Pbcrit<-Pb[c(ceiling(ll*nboot),ceiling(sl*nboot)),,]

for ( i in 1 : length(svars)) {
    svars[[i]]<-which(sign(drop(Pbcrit[1,i]*Pbcrit[2,i]))==1)
 } #

return(list(Pb, Pbcrit,svars))
}
#
# ############################################################
orth.proc <-function(x,y) {
#computes orthogonal procrustes rotation projecting matrix y onto the
#subspace spanned by matrix x.
# syntax: [r,yrot]<-orth.proc(x,y)
# where r is the rotation matrix and yrot is the procrustes rotated y

  rsvd<-svd(t(y)*x, 0)

r<-rsvd$u*t(rsvd$v)
yrot<-y*r

return(list(r,yrot))
 }

