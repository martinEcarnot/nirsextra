asca <-function(...)  {

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
} elseif (nargin==2) {
    X<-varargin{1}
    desmatrix<-varargin{2}
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
} elseif (nargin==3) {
    X<-varargin{1}
    desmatrix<-varargin{2}
    opt<-varargin{3}
 } #


#preprocessing
ntot<-nrow(X)

switch(opt$preproc,
    none={Xp<-X},
    auto={Xp<-X./repmat(std(X), ntot, 1)},
    pareto={Xp<-X./repmat(sqrt(std(X)),ntot, 1)},
)

model$Xdata$PreprData<-Xp
model$Xdata$TotalSSQ<-sum(sum(Xp^2))
model$Design<-desmatrix

dmain<-gendmat(desmatrix)
#[dmatrices,desterms,deslabels, desorder]<-createdesign(dmain)
ldpatrices = createdesign(dmain)

Xd<-Xp

for ( i in 1 : length(dmatrices)) {
    Xeff<-dmatrices{i}*pinv(dmatrices{i})*Xd
    ssqEff<-sum(sum(Xeff.^2))
    Xd<-Xd-Xeff
    l<-strtrim(deslabels{i})
    eval(paste0('model.X', l,'.EffectMatrix<-Xeff'))
    eval(paste0('model.X', l,'.EffectSSQ<-ssqEff'))

    if (i==1) {
        ssqtot<-sum(sum(Xd.^2))
        model.Xdata.CenteredData<-Xd
        model.Xdata.CenteredSSQ<-ssqtot
    } else {
        expVar<-100*(ssqEff/ssqtot)
        eval(paste0('model.X', l,'.EffectExplVar<-expVar'))
     } #

    eval(paste0('model.X', l,'.DesignMatrix<-dmatrices{i}'))
    eval(paste0('model.X', l,'.DesignTerms<-desterms{i}'))
    eval(paste0('model.X', l,'.EffectLabel<-strtrim(deslabels{i})'))
    eval(paste0('model.X', l,'.TermOrder<-desorder(i)'))

 } #

model$XRes$EffectMatrix<-Xd
model$XRes$EffectSSQ<-sum(sum(Xd.^2))
model$XRes$EffectExplVar<-100*(model.XRes.EffectSSQ/ssqtot)
model$XRes$EffectLabel<-'Res'
model$XRes$TermOrder<-max(desorder)+1

for ( i in 2 : length(dmatrices)) {
    l<-strtrim(deslabels{i})

    if (ischar(opt.reducedmodel) && strcmp(opt.reducedmodel, 'standard')) {
        remfact<-desterms(c(2:i-1,i+1:end))
    } else {
        remfact<-opt.reducedmodel{i}
     } #
    Xx<-model.Xdata.CenteredData
    for ( j in 1 : length(remfact)) {
        m<-strtrim(char(64+remfact{j}))
        eval(paste0('Xx<-Xx-model.X',m,'.EffectMatrix'))
     } #
    eval(paste0('model.X', l,'.ReducedMatrix<-Xx'))
 } #

model.TermLabels<-deslabels
model.Options<-opt


#Permutation tests, if required
model<-ascaptest(model)

#SCA Models
model<-ascasca(model)

#bootstrap, if required
model<-ascaboot(model)


}






#################################

createdesign <-function(designmat) {

ns<-nrow(designmat)
nfact<-ncol(designmat)

nfact<-length(designmat)
ns<-nrow(designmat{1})

indmat<-fullfact(repmat(2,1,nfact))
nmat<-size(indmat,1)
dmatrices<-cell(1,nmat)
desterms<-cell(1,nmat)
deslabels<-cell(1,nmat)
desorder<-matrix(0,nmat,1)



for ( i in 1 : nmat) {
    dm<-matrix(1,ns,1)
    for ( j in 1 : nfact) {
        if (indmat(i,j)==1) {
            effmat<-designmat{j}
        }else{
            effmat<-matrix(1,ns,1)
         } #
        dm<-kron(dm,effmat)
        dm<-dm(1:ns+1:end) #,:)
     } #
    dmatrices{i}<-dm
    desterms{i}<-find(indmat(i,)==1)
    deslabels{i}<-char(64+find(indmat(i,)==1))
    desorder(i)<-length(desterms{i})


 } #

deslabels<-sort(char(deslabels),2)

#Sorting according to increasing order of interactions
#[deslabels, newindex]<-sortrows(deslabels)
ldeslabels = sortrows(deslabels)

desterms<-desterms(newindex)
dmatrices<-dmatrices(newindex)
desorder<-desorder(newindex)
deslabels<-t(cellstr(deslabels))
# deslabels{1}<-'Mean'
deslabels[1]<-'Mean'

return(list(dmatrices,desterms, deslabels, desorder))
}


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

#############################################
ascaptest <-function(ascamodel) {
pmodel

pmodel<-ascamodel
dlab<-ascamodel.TermLabels

if strcmp(ascamodel.Options.permtest, 'on') {

    signfacts<-cell(length(dlab)-1,1)
    sc<-0
}else{
    signfacts<-[]
 } #


for ( i in 2 : length(dlab)) {
    l<-strtrim(dlab{i})
    if strcmp(ascamodel.Options.permtest, 'on') {
        if (ischar(ascamodel.Options.permfacts) && strcmp(ascamodel.Options.permfacts, 'all')) {

            eval(['Xr<-ascamodel.X', l, '.ReducedMatrix '])
            eval(['Dr<-ascamodel.X', l, '.DesignMatrix '])
            ssqp<-ptest(Xr,Dr, ascamodel.Options.nperm)
            eval(['seff<-ascamodel.X', l, '.EffectSSQ '])
            p<-length(find(ssqp><-seff))./ascamodel.Options.nperm
            if (p<<-0.05) {
                sc<-sc+1
                signfacts{sc}<-l
             } #


            eval(['pmodel.X',l,'.EffectSignif.NullDistr<-ssqp'])
            eval(['pmodel.X',l,'.EffectSignif.p<-p'])
        }else{
            if (ismember(l, ascamodel.Options.permfacts)) {
                eval(['Xr<-ascamodel.X', l, '.ReducedMatrix '])
                eval(['Dr<-ascamodel.X', l, '.DesignMatrix '])
                ssqp<-ptest(Xr,Dr, ascamodel.Options.nperm)
                eval(['seff<-ascamodel.X', l, '.EffectSSQ '])
                p<-length(find(ssqp><-seff))./ascamodel.Options.nperm
                if (p<<-0.05) {
                    sc<-sc+1
                    signfacts{sc}<-l
                 } #

                eval(['pmodel.X',l,'.EffectSignif.NullDistr<-ssqp'])
                eval(['pmodel.X',l,'.EffectSignif.p<-p'])
            }else{
                eval(['pmodel.X',l,'.EffectSignif.NullDistr<-[]'])
                eval(['pmodel.X',l,'.EffectSignif.p<-[]'])


             } #
         } #
    }else{
        eval(['pmodel.X',l,'.EffectSignif.NullDistr<-[]'])
        eval(['pmodel.X',l,'.EffectSignif.p<-[]'])
     } #
 } #

if (strcmp(ascamodel.Options.permtest, 'on')) {
    signfacts<-signfacts(1:sc)
 } #

pmodel.SignificantTerms<-signfacts

}

########################################################
ptest <-function(X,D, nperm) {

ns<-size(X,1)  #Number of samples
ssqp<-matrix(0,nperm,1)   #Initialization of the permuted SSQ vector

for ( i in 1 : nperm) {
    hh<-randperm(ns)
    Xpp<-D(hh,:)*pinv(D(hh,:))*X
    ssqp(i)<-sum(sum(Xpp.^2))
 } #

return(ssqp)
}

##########################################################
ascasca <-function(ascamodel) {

smodel<-ascamodel
dlab<-ascamodel.TermLabels

for ( i in 2 : length(dlab)) {
    l<-strtrim(dlab{i})
    eval(['Xr<-ascamodel.X', l, '.EffectMatrix '])
    eval(['ssqr<-ascamodel.X', l, '.EffectSSQ '])

    R<-rank(Xr)
    [u,s,P]<-svds(Xr,R)
    t<-u*s
    taug<-(Xr+ascamodel.XRes.EffectMatrix)*P
    varex<-100*(diag(s).^2)/ssqr

    eval(['smodel.X', l, '.SCA.Model.Scores<-t '])
    eval(['smodel.X', l, '.SCA.Model.ScoreswithRes<-taug '])
    eval(['smodel.X', l, '.SCA.Model.Loadings<-P '])
    eval(['smodel.X', l, '.SCA.Model.ExplVar<-varex '])
#       eval(['smodel.X', l, '.SCA.Model.s<-s ']) # Only if we want to
#       trace error : Q residuals etc...
 } #

return(smodel)

}

##########################################################
ascaboot <-function(ascamodel) {

bmodel<-ascamodel
dlab<-ascamodel.TermLabels
Xd<-ascamodel.Xdata.CenteredData

for ( i in 2 : length(dlab)) {
    l<-strtrim(dlab{i})

    if (strcmp(ascamodel.Options.bootstrap, 'all')) {
        if (strcmp(ascamodel.Options.bootmatrix, 'original')) {

            Xd<-ascamodel.Xdata.CenteredData
        }elseif (strcmp(ascamodel.Options.bootmatrix, 'reduced')) {
            eval(['Xd<-ascamodel.X', l, '.ReducedMatrix '])
         } #

        eval(['Dd<-ascamodel.X', l, '.DesignMatrix '])
        eval(['Pd<-ascamodel.X', l, '.SCA.Model.Loadings '])
        [Pb, Pbcrit, svars]<-bootload(Xd,Dd, Pd, ascamodel.Options.confl, ascamodel.Options.nboot)
    } elseif (strcmp(ascamodel.Options.bootstrap, 'signif')) {
        if (ismember(l, ascamodel.SignificantTerms)) {
            if (strcmp(ascamodel.Options.bootmatrix, 'original')) {

                Xd<-ascamodel.Xdata.CenteredData
            } elseif (strcmp(ascamodel.Options.bootmatrix, 'reduced')) {

                eval(['Xd<-ascamodel.X', l, '.ReducedMatrix '])
             } #
            eval(['Dd<-ascamodel.X', l, '.DesignMatrix '])
            eval(['Pd<-ascamodel.X', l, '.SCA.Model.Loadings '])
            [Pb, Pbcrit, svars]<-bootload(Xd,Dd, Pd, ascamodel.Options.confl, ascamodel.Options.nboot)
        } else {
            Pb<-[]
            Pbcrit<-[]
            svars<-[]
         } #


    } elseif (strcmp(ascamodel.Options.bootstrap, 'off')) {
        Pb<-[]
        Pbcrit<-[]
        svars<-[]
     } #


    switch (ascamodel.Options.bootsave,
        'all'={
            eval(['bmodel.X', l, '.SCA.Bootstrap.Loadings<-Pb '])
            eval(['bmodel.X', l, '.SCA.Bootstrap.ConfIntervals<-Pbcrit '])
            eval(['bmodel.X', l, '.SCA.Bootstrap.SignificantVariables<-svars '])
}
        'confint'={
            eval(['bmodel.X', l, '.SCA.Bootstrap.Loadings<-[] '])
            eval(['bmodel.X', l, '.SCA.Bootstrap.ConfIntervals<-Pbcrit '])
            eval(['bmodel.X', l, '.SCA.Bootstrap.SignificantVariables<-svars '])
}
        'signvars'={
            eval(['bmodel.X', l, '.SCA.Bootstrap.Loadings<-[] '])
            eval(['bmodel.X', l, '.SCA.Bootstrap.ConfIntervals<-[] '])
            eval(['bmodel.X', l, '.SCA.Bootstrap.SignificantVariables<-svars '])
}
     ) #

 } #

}


#########################################################
bootload <-function(Xd,Dd, Pd, confl, nboot) {
#Bootstrap

sl<-(confl+1)/2
ll<-1-sl
svars<-cell(1,size(Pd,2))

Pb<-matrix(0,[nboot size(Pd)])

bootp<-matrix(0,size(Xd,1),1)
lev<-unique(Dd, 'rows')

for ( i in 1 : nboot) {
    for ( j in 1 : nrow(lev,1)) {
        xx<-find(ismember(Dd, lev(j,:), 'rows')==1)
        yy<-ceil(length(xx)*rand(1,length(xx)))
        bootp(xx)<-xx(yy)
     } #
    Xboot<-Xd(bootp,:)
    Xdp<-Dd*pinv(Dd)*Xboot
    [~,~,vb]<-svds(Xdp,size(Pd,2))
    [~,Pb(i,:,:)]<-orth.proc(Pd,vb)
 } #
Pb<-sort(Pb)
Pbcrit<-Pb([ceil(ll*nboot) ceil(sl*nboot)],:,:)

for ( i in 1 : length(svars)) {
    svars{i}<-find(sign(squeeze(Pbcrit(1,:,i).*Pbcrit(2,:,i)))==1)
 } #

return(list([Pb, Pbcrit,svars]))
}

############################################################
orth.proc <-function(x,y) {
#computes orthogonal procrustes rotation projecting matrix y onto the
#subspace spanned by matrix x.
# syntax: [r,yrot]<-orth.proc(x,y)
# where r is the rotation matrix and yrot is the procrustes rotated y

[u,~,v]<-svd(y'*x, 0)
r<-u*v'
yrot<-y*r

return(list(r,yrot)))
 }

