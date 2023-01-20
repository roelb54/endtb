age.group.from.lower.age.limits = function(age.from,lower.age.limits){
  if (is.null(lower.age.limits)){
    age.group = rep("ALL",length(age.from))
  }else{
    lower=lower.age.limits
    upper=c(lower.age.limits[-1],200)
    age.group = rep(NA,length(age.from))
    for (i in seq_along(lower)){
      age.group[age.from %in% lower[i]:upper[i]]=paste0("A",lower[i])
    }
  }
  age.group
}
validate.target.states=function(target=NULL,names=NULL,p=NULL){
  names(target)=names
  substrs = stri_split_fixed(target,",")
  for (i in seq_along(names)){
    check = substrs[[i]] %in% p$DIMNAMESLIST[[names(target)[i]]]
    assert_that(all(check),msg = paste("target states",paste(substrs[[i]][!check], collapse=","),"not found in",paste(p$DIMNAMESLIST[[names(target)[i]]],collapse=",")))
  }
}
query.from.target=function(target=NA,p=NULL){ # NOTE the dead are excluded in the function calling this fn
  n = length(target)
  assert_that(n==8,msg="target row should be 8 elements")
  target[] = stri_replace_all_regex(target,replacement = "",pattern="[:space:]")
  assert_that(!is.na(target[1]),msg="NA not allowed in year column")
  mult = rep(F,n)
  mult[c(2:6)]=T # i.e. multiple states
  num = rep(F,n)
  num[c(1,7,8)]=T # i.e. year, age_from, age_thru
  sel = mult & (target=="all" | target=="NA" | is.na(target))
  target[sel]=NA
  NAs = is.na(target)
  target[NAs]="" 
  validate.target.states(target[!NAs & mult],names(target)[!NAs & mult],p)
  target[!NAs & mult] = stri_replace(target[!NAs & mult],"','",regex=",",mode="all")
  target[!NAs & mult] = stri_replace(target[!NAs & mult],"c('",regex="^")
  target[!NAs & mult] = stri_replace(target[!NAs & mult],"')",regex="$")
  target[!NAs & mult] = paste0(" %in% ",target[!NAs & mult])
  target[!NAs & mult] = paste(names(target)[!NAs & mult],target[!NAs & mult])
  target[1] = paste0(names(target)[1],"==",target[1])
  if (target[7]==0 & target[8]==99){
    target[c(7,8)] = paste0(names(target)[c(7,8)],"==",target[c(7,8)])
  }else{
    target[7] = paste0(names(target)[7],">=",target[7])
    target[8] = paste0(names(target)[8],"<=",target[8])
  }
  s = stri_join(target[target!=""],collapse=" & ")
  # print(s)
  return(s)
}

read.targets = function(fname){
  df = read.csv(fname,header=T,stringsAsFactors = F, fileEncoding = 'UTF-8-BOM')
  assert_that(length(unique(df[,1]))==1,msg="a target file should contain only one country")
  df = df[,-1]
  #df = as.data.frame(sapply(df, function (x) gsub ("\"", "", x)))
  #df = as.data.frame(sapply(df, function (x) gsub ("\'", "", x)))
  assert_that(sum(is.na(df$value))==0, msg="NAs not allowed in value column of targets file")
  assert_that(sum(is.na(df$lo)   )==0, msg="NAs not allowed in lo column of targets file")
  assert_that(sum(is.na(df$hi)   )==0, msg="NAs not allowed in hi column of targets file")
  assert_that(sum(is.na(df$err)  )==0, msg="NAs not allowed in err column of targets file")
  assert_that(all((df$hi-df$lo)>0), msg="not all hi > lo")
  assert_that(all(df$value>=df$lo), msg="not all values >= lo")
  assert_that(all(df$value<=df$hi), msg="not all values <= hi")
  assert_that(all(df$err>0 | df$err==-1), msg="not all err > 0 or -1")
  df  
}

exclude.dead=function(z){
  hivdead  = grepl("dead$",z$HIV,ignore.case = T) 
  tbdead   = grepl("dead$",z$TB,ignore.case = T) 
  riskdead = grepl("dead$",z$RISK,ignore.case = T) 
  z[!(hivdead | tbdead | riskdead),]
}
calc.stock=function(z,target,p=NULL){
  assert_that(target$type=="stock",msg="stock type should be 'stock' and stock cause should be NA")
  q = query.from.target(target[c(1:6,9,10)],p)
  modlog(level='DEBUG',msg=paste0("querying stock in target evaluation with: ",q))
  as.numeric(z[eval(parse(text=q)),.(value=sum(value))])
}

calc.flow=function(z,target,p=NULL){
  causes = stri_trim_both(unlist(strsplit(target$cause,",")))
  allowed_causes = c('Tm',paste0(names(p$inci.dim.names[!is.null(p$inci.dim.names)]),'Xi'),paste0(p$DIMNAMES,'p'),paste0('TBtr_',names(p$TBtr)),paste0('HIVtr_',names(p$HIVtr)))
  err = sum(causes %in% allowed_causes)!=length(causes)
  assert_that(!err, msg="undefined cause for flow")
  if (err){
    modlog(level="ERROR",msg="undefined cause for flow")
  }
  if (target$type=="inflow"){ inflow = T }else if (target$type=="outflow"){ inflow = F }else{ inflow = NULL }
  assert_that(!is.null(inflow), msg="type for flow should be 'inflow' or 'outflow'")
  q = query.from.target(target[c(1:6,9,10)],p)
  q = paste0(q," & flow ==",ifelse(inflow,"'in'","'out'"))
  q = paste0(q," & subject %in% c('",paste0(causes,collapse="','"),"')")
  modlog(level='DEBUG',msg=paste0("querying flow in target evaluation with: ",q))
  as.numeric(z[eval(parse(text=q)),.(value=sum(value))])
}

calc.targets = function(z,targets,p=NULL){
  model = rep(NA, nrow(targets))
  for (i in seq_along(targets$year_from)){
    x_range = c(1, 3:9, 21,22)    ; x_row   = targets[i,x_range] ; names(x_row)[1]="year" ; names(x_row)[2:8]=substring(names(x_row)[2:8],2)
    y_range = c(1,13:19,21,22)    ; y_row   = targets[i,y_range] ; names(y_row)[1]="year" ; names(y_row)[2:8]=substring(names(y_row)[2:8],2)
    thruna = is.na(targets$year_thru) ; targets$year_thru[thruna]=targets$year_from[thruna]
    model[i]=0
    if (targets$year_thru[i]>targets$year_from[i]){ # i.e. a range of years
      for (yr in targets$year_from[i]:targets$year_thru[i]){
        x_row$year = y_row$year = yr
        if     (x_row$type == "stock")  {x = calc.stock(exclude.dead(z$stocks),x_row,p=p)}
        else if(x_row$type == "inflow") {x = calc.flow(z$flows,x_row,p=p)}
        else if(x_row$type == "outflow"){x = calc.flow(z$flows,x_row,p=p)}
        if (is.na(targets$a[i])){ a = 1} else {a = targets$a[i]}
        if (is.na(targets$b[i])){ b = 0} else {b = targets$b[i]}
        if (is.na(targets$c[i]) | yr==targets$year_from[i])  { c = 0 } else {c = targets$c[i]}
        if     (y_row$type == "stock")  {y = calc.stock(exclude.dead(z$stocks),y_row,p=p)}
        else if(y_row$type == "inflow") {y = calc.flow(z$flows,y_row,p=p)}
        else if(y_row$type == "outflow"){y = calc.flow(z$flows,y_row,p=p)}
        eval(parse(text=paste0("model[i]=model[i]+",targets[i,]$fn)))
      }
    }else{
      if (is.na(targets$a[i])){ a = 1} else {a = targets$a[i]}
      if (is.na(targets$b[i])){ b = 0} else {b = targets$b[i]}
      if (is.na(targets$c[i])){ c = 0} else {c = targets$c[i]}
      if     (x_row$type == "stock")  {x = calc.stock(exclude.dead(z$stocks),x_row,p=p)}
      else if(x_row$type == "inflow") {x = calc.flow(z$flows,x_row,p=p)}
      else if(x_row$type == "outflow"){x = calc.flow(z$flows,x_row,p=p)}
      if     (y_row$type == "stock")  {y = calc.stock(exclude.dead(z$stocks),y_row,p=p)}
      else if(y_row$type == "inflow") {y = calc.flow(z$flows,y_row,p=p)}
      else if(y_row$type == "outflow"){y = calc.flow(z$flows,y_row,p=p)}
      eval(parse(text=paste0("model[i]=",targets[i,]$fn)))
    }
  }
  # subset = targets$err>=0
  # weighted_residuals = residuals = fit = rep(NA,length(targets$err))
  # fit[subset] = model[subset] >= targets[subset,]$lo & model[subset] <= targets[subset,]$hi
  # residuals[subset] = model[subset] - targets[subset,]$value
  # weighted_residuals[subset] = residuals[subset] / targets[subset,]$err
  cbind(targets,model=model)
}

eval.output.vs.targets = function(output,targets,p=NULL,exclude.intermediates=F){
  
  okfrom  = all(targets$year_from %in% unique(output$stocks$year))
  selthru = !is.na(targets$year_thru)
  okthru  = all(targets$year_thru[selthru] %in% unique(output$stocks$year))
  fromyrs = targets$year_from[selthru]
  thruyrs = targets$year_thru[selthru]
  
  okranges = T
  if (sum(selthru>0)){
    for (i in seq_along(fromyrs)){
      yrs = fromyrs[i]:thruyrs[i]
      okranges = okranges & all(yrs %in% unique(output$stocks$year))
    }
  }
  ok = okfrom & okthru & okranges
  if (!ok){
    if (!okfrom){
      cat("target year_from specified that are not in output:\n")
      indxs = !(targets$year_from %in% unique(output$stocks$year))
      apply(targets[indxs,],1,function(x) cat("year_from:",x["year_from"]," name:",x["name"],"\n"))
    }
    if (!okthru){
      cat("target year_thru specified that are not in output:\n")
      indxs = !(targets$year_thru %in% unique(output$stocks$year))
      apply(targets[indxs,],1,function(x) cat("year_thru:",x["year_thru"]," name:",x["name"],"\n"))
    }
    if (!okranges){
      cat("target year range specified that are not in output:\n")
      for (i in seq_along(fromyrs)){
        yrs = fromyrs[i]:thruyrs[i]
        if (!all(yrs %in% unique(output$stocks$year))){
          cat("not all years of range :",yrs," in stocks output \n")
        }
      }
    }
  }
  assert_that(ok,msg="targets specified that are not in output, modify xml file")
  selx = (targets$xtype %in% c("stock","inflow","outflow"))
  sely = (targets$ytype %in% c("stock","inflow","outflow"))
  assert_that(all(selx == sely),msg="inconsistent definition of combined targets")
  
  result = calc.targets(output,targets[selx,],p=p)
  
  result2 = NULL
  subtargets = targets[!selx,]
  if (nrow(subtargets)>0){
    n = nrow(subtargets)
    for (i in 1:n){
      xtype = ytype = xcause = ycause = NA
      row  =  subtargets[i,]
      xtype = result[  result$name == row$xtype & result$year_from == row$year_from,"model"]
      ytype = result[  result$name == row$ytype & result$year_from == row$year_from,"model"]
      if (!(is.na(row$xcause) | row$xcause=="")){
        xcause = result[  result$name == row$xcause & result$year_from == row$year_from,"model"]
      }
      if (!(is.na(row$ycause) | row$ycause=="")){
        ycause = result[  result$name == row$ycause & result$year_from == row$year_from,"model"]
      }
      eval(parse(text=paste0("row$model=",row$fn)))
      result2 = rbind(result2,row)
    }
  }
  z = rbind(result[result$err>=0,], result2)
  weighted_residuals = residuals = fit = z$err * 0
  fit = z$model >= z$lo & z$model <= z$hi
  residuals = z$model - z$value
  weighted_residuals = residuals / z$err
  z = cbind(z,fit,residuals,weighted_residuals)
  modlog(level="DEBUG",msg="successful call of eval.output.vs.targets()")
  z
#  if (exclude.intermediates){
#    return(rbind(result[!is.na(result$fit),],result2))
#  }else{
#    return(rbind(result,result2))
#  }
}

