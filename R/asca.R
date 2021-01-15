asca <-function(...)  {

library(MASS)

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
  # #############################################
  ascaptest <-function(ascamodel) {

    pmodel<-ascamodel
    dlab<-ascamodel$TermLabels

    if  (ascamodel$Options$permtest=='on') {

      signfacts<-list(length(dlab)-1)
      sc<-0
    }else{
      signfacts<-NULL
    }

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
          }

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
            }

            eval(parse(text=paste0('pmodel$X',l,'$EffectSignif$NullDistr<-ssqp')))
            eval(parse(text=paste0('pmodel$X',l,'$EffectSignif$p<-p')))
          }else{
            eval(parse(text=paste0('pmodel$X',l,'$EffectSignif$NullDistr<-NULL]')))
            eval(parse(text=paste0('pmodel$X',l,'$EffectSignif$p<-NULL')))
          }
        }
      }else{
        eval(parse(text=paste0('pmodel$X',l,'$EffectSignif$NullDistr<-[]')))
        eval(parse(text=paste0('pmodel$X',l,'$EffectSignif$p<-[]')))
      }
    }

    if (ascamodel$Options$permtest == 'on') {
      signfacts<-signfacts[1:sc]
    }

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
      t<-rsvd$u %*% rsvd$d[1:R]
      taug<-(Xr+ascamodel$XRes$EffectMatrix) %*% rsvd$v[,1:R]
      varex<-100*(diag(as.matrix(rsvd$d[1:R]))^2)/ssqr

      eval(parse(text=paste0('smodel$X', l, '$SCA$Model$Scores<-t ')))
      eval(parse(text=paste0('smodel$X', l, '$SCA$Model$ScoreswithRes<-taug ')))
      eval(parse(text=paste0('smodel$X', l, '$SCA$Model$Loadings<-rsvd$v[,1:R] ')))
      eval(parse(text=paste0('smodel$X', l, '$SCA$Model$ExplVar<-varex ')))
      #       eval(parse(text=paste0('smodel$X', l, '$SCA$Model$s<-s ']) # Only if we want to
      #       trace error : Q residuals etc...
    }
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
        }

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
        }

      } else if (ascamodel$Options$bootstrap == 'off') {
        Pb<-NULL
        Pbcrit<-NULL
        svars<-NULL
      }

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
                eval(parse(text=paste0('bmodel$X', l, '$SCA$Bootstrap$Loadings<-NULL ')))
                eval(parse(text=paste0('bmodel$X', l, '$SCA$Bootstrap$ConfIntervals<-Pbcrit ')))
                eval(parse(text=paste0('bmodel$X', l, '$SCA$Bootstrap$SignificantVariables<-svars ')))
              },
              signvars={
                eval(parse(text=paste0('bmodel$X', l, '$SCA$Bootstrap$Loadings<-NULL ')))
                eval(parse(text=paste0('bmodel$X', l, '$SCA$Bootstrap$ConfIntervals<-NULL ')))
                eval(parse(text=paste0('bmodel$X', l, '$SCA$Bootstrap$SignificantVariables<-svars ')))
              }
      )
    }
    return(bmodel)
  }
  # #########################################################
  bootload <-function(Xd,Dd, Pd, confl, nboot) {
    #Bootstrap

    sl<-(confl+1)/2
    ll<-1-sl
    if (is.null(dim(Pd))) {dim(Pd)=c(length(Pd),1)}
    if (is.null(dim(Dd))) {dim(Dd)=c(length(Dd),1)}
    svars<-vector("list", dim(Pd)[2])

    Pb<- array(0, c(nboot, dim(Pd)[1], dim(Pd)[2])) # Pb=zeros([nboot size(Pd)]);

    bootp<-matrix(0,nrow(Xd))
    lev<-unique(Dd) #unique(Dd, 'rows')

    for ( i in 1 : nboot) {
      for ( j in 1 : nrow(lev)) {
        xx<-which(!is.na(match(data.frame(t(Dd)), data.frame(lev[j,]))))
        yy<-ceiling(length(xx)*runif(length(xx)))  # rand(1,length(xx)))
        bootp[xx]<-xx[yy]
      } #
      Xboot<-Xd[bootp,]
      Xdp<-Dd %*% ginv(Dd) %*% Xboot
      rsvd<-svd(Xdp,dim(Pd)[2])  # [~,~,vb]<-svds(Xdp,size(Pd,2))
      resOrth.proc<-orth.proc(Pd,rsvd$v[,1:dim(Pd)[2]])
      Pb[i,,]=resOrth.proc$yrot
    } #
    ## Pb<-sort(Pb)
    # Pbdim0=dim(Pb)
    # dim(Pb)=c(Pbdim0[1],Pbdim0[2]*Pbdim0[3])
    # Pb=apply(Pb,2,sort)
    # dim(Pb)=Pbdim0

    Pbcrit<-Pb[c(ceiling(ll*nboot),ceiling(sl*nboot)),,]

    for ( i in 1 : length(svars)) {
      svars[[i]]<-which(sign(drop(Pbcrit[1,,i]*Pbcrit[2,,i]))==1)
    } #

    return(list(Pb=Pb, Pbcrit=Pbcrit,svars=svars))
  }
  # ############################################################
  orth.proc <-function(x,y) {
    #computes orthogonal procrustes rotation projecting matrix y onto the
    #subspace spanned by matrix x.
    # syntax: [r,yrot]<-orth.proc(x,y)
    # where r is the rotation matrix and yrot is the procrustes rotated y

    rsvd<-svd(t(y) %*% x)

    r<-rsvd$u*t(rsvd$v)
    yrot<-y %*% r

    return(list(r=r,yrot=yrot))
  }
  ############################"

  args = list(...)

if (nargs()==1 && args[[1]]=='options') {
  # Extract default options
    opt=list()
    opt$preproc<-'none'
    opt$reducedmodel<-'standard'
    opt$permtest<-'on'
    opt$permfacts<-'all'
    opt$nperm<-100
    opt$bootstrap<-'off'
    opt$bootmatrix<-'original'
    opt$bootsave<-'confint'
    opt$nboot<-100
    opt$plots<-'on'
    opt$confl<-0.95
    model<-opt
    return(model)
} else if (nargs()==2) {
    X<-args[[1]]
    desmatrix<-args[[2]]
    opt=list()
    opt$preproc<-'none'
    opt$reducedmodel<-'standard'
    opt$permtest<-'on'
    opt$permfacts<-'all'
    opt$nperm<-100
    opt$bootstrap<-'off'
    opt$bootmatrix<-'original'
    opt$bootsave<-'confint'
    opt$nboot<-100
    opt$plots<-'on'
    opt$confl<-0.95
} else if (nargs()==3) {
    X<-args[[1]]
    desmatrix<-args[[2]]
    opt<-args[[3]]
 }


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
     }
    eval(parse(text=paste0('model$X', l,'$DesignMatrix<-dmatrices[[i]]')))
    eval(parse(text=paste0('model$X', l,'$DesignTerms<-desterms[[i]]')))
    eval(parse(text=paste0('model$X', l,'$EffectLabel<-deslabels[[i]]')))
    eval(parse(text=paste0('model$X', l,'$TermOrder<-desorder[i]')))
 }

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
     }
    Xx<-model$Xdata$CenteredData
    for ( j in 1 : length(remfact)) {
        m<-paste0(LETTERS[remfact[[j]]],collapse ="")
        eval(parse(text=paste0('Xx<-Xx-model$X',m,'$EffectMatrix')))
     }
    eval(parse(text=paste0('model$X', l,'$ReducedMatrix<-Xx')))
 }

model$TermLabels<-deslabels
model$Options<-opt


#Permutation tests, if required
model<-ascaptest(model)

#SCA Models
model<-ascasca(model)

#bootstrap, if required
model<-ascaboot(model)

}
