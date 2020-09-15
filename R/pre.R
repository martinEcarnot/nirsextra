pre <- function(sp,p)
{

  # Pretreatment
  # sp =spectra, p= list of pretreatment

  # Liste de prtraitements	Description			Arguments
  #	'snv'		Sdandard Normal Variate			''
  # 'msc'		Multi-scater Correction			ref :1 spectre  partir duquel est calcul la correction (dfaut= mean(Xi.d))
  #	'ma'		Moving average			wsz: taille de la fenetre glissante (dfaut=3)
  #	'red'		Reduction du nbre de l.o.			[d,gap,f]: nombre de longueurs d'onde passes au dbut, pas d'chantillonnage, nombre de longueurs d'onde passes  la fin
  #	'der'		derivee				ordre de drivation  (dfaut=1)
  #	'sder'		derivee SA1SIR			[polynom_order,window_size,derivative_order]   (dfaut=[3,5,1])
  #	'adj'		ajustement des sauts de detecteurs		Longueurs d'onde auxquelles la correction doit etre faite (defaut=[651,1451])
  #	'detr'		detrend				Degr du polynome d'ajustement pour le detrend
  # 'vs'		Variance Scaling			ref :ecart-type du C-Set de reference (dfaut= std(Xi.d))
  #	'norm_rms'	Normalalisation par la moy quadratique	''
  #	'dec'		decalage des l.o.			intervale dans lequel le minimum sera trouve (def=[501 576])
  #	'ref2abs'	Passe de reflectance à l'absorbance	''
  # 'offset' rajoute un offset (def=0)
  # 'v01' Cale les spectres entre 0 et 1
  # 'del'   Enleve des morceaux de spectres (pas seulement aux bords)

  #Effectue les prtraitements sur les spectres

  # browser()
  if (missing(p)) {p=list()}
  if (length(p)<3) {dim(p)=c(1,2)}
  nt=nrow(p)
  pf=list()

  for (i in 1:nt) {
  str='Rien'

  if (p[[i,1]] == 'snv') {
    sp = t(scale(t(sp)))
    str = ''
  }

  if (p[[i,1]] == 'vs') {
    if (p[[i,2]] == '') {
      ref = apply(sp, 2, sd)
      str = 'std'
    }
    else {
      ref = p[[i,2]]
      str = 'udef'
    }
    sp = sp / rep(1, nrow(sp)) %*% t(ref)
  }

  if (p[[i,1]] == 'msc') {
    if (p[[i,2]] == '') {
      sp = msc(sp)
      str = mscmsc(matrix(rep(NA, ncol(sp)))) # create empty msc object
      mr_attr = attributes(mr)
      dim(mr_attr) = NULL
      attributes(str) = mr_attr
    } # fill empty msc model with good attributes (saves memory)
    else
      sp = predict(str, sp)
    str = p[[i,2]]
  }


  if (p[[i,1]] == 'red') {
    if (length(p[[i,2]]) == 1) {
      str = c(0, 0, 1)
    }
    else{
      str = p[[i,2]]
      sp = sp[, seq(str[1], ncol(sp) - str[2], str[3])]
    }
  }

  if (p[[i,1]] == 'del') {
    if (length(p[[i,2]]) == 1) {
      str = c(0, 0)
    }
    else{
      str = p[[i,2]]
      sp = sp[, -(str[1]:str[2])]
    }
  }

  if (p[[i,1]] == 'sder') {
    if (length(p[[i,2]]) == 1) {
      str = c(0, 1, 1)
    }
    else{
      str = p[[i,2]]
      sp = savgol(sp, m = str[1], p = str[2], n = str[3])
    }
  }



  if (p[[i,1]] == 'adj') {
    if (length(p[[i,2]]) == 1) {
      str = c(651, 1451)
    }
    else{
      str = p[[i,2]]
    }
      sp = adj_asd(sp, str)

  }


  if(p[[i,1]]=='detr') {
    str = p[[i,2]]
    eval(parse(text=paste0("sp=detrend(sp,",str,")")))
  }

  if(p[[i,1]]=='ref2abs') {
    sp=-log(sp)
    str = p[[i,2]]
  }

  #
  #   if(p[[i,1]]=='norm_rms') {
  #   Xi=norm_rms(Xi);
  #   str='';
  #
  #   if(p[[i,1]]=='dec') {
  #   if strcmp(p{i,2},'') == 1; interval=[501 576]; else interval=p{i,2}; end
  #   Xi=rec_lo(Xi,interval);
  #   str='';
  #

  #
  #   if(p[[i,1]]=='abs2ref') {
  #   Xi.d=abs2ref(Xi.d);
  #   str='';
  #
  #   if(p[[i,1]]=='log'
  #   Xi.d=log_mec(Xi.d);
  #   str='';
  #
  #   if(p[[i,1]]=='offset') {
  #   if strcmp(p{i,2},'') == 1; offs=0; else offs=p{i,2}; end
  #   Xi.d=Xi.d+offs;
  #   str=num2str(offs);
  #
  #   if(p[[i,1]]=='v01') {
  #   mmax=max(Xi.d');
  # Xi.d=Xi.d./repmat(mmax',1,size(Xi.d,2));
  #                   str='';
  #
  #                   end#switch
  #
  #                   pf=[pf, p{i,1}, '-',str,'-'];
  #
  #                   end
  #
  #                   if nt == 0;
  #                   pf='Rien';
  #                   end
  #                   Xf=Xi;
pf=rbind(pf,list(p[[i,1]],str))
  }
  attr(sp,"pre")=pf
return(sp)
}
