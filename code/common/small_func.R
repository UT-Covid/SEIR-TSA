library(lubridate)
da = function(x, year=2020){
        m = as.integer( floor(x) )
        d = as.integer(round( (x-m)*100 ))
        ymd( sprintf("%d-%02d-%02d",year,m,d))
}

daS = function( from, to=c(), by=c(),len=c(), ...) {
    x = c(from, to)
    x=sort(x)
    if( is.numeric(x) ) {
        if( x[1]<0) {
            dx = c( da(x[2])+x[1], da(x[2]) )
        } else if( length(x)>1  & x[2] == floor(x[2]) ) {
            dx = c( da(x[1]), da(x[1])+x[2] )
        } else {
            dx = da(x)
        }
    } else {
        dx = x
    }
    if( !is.null(by) & !is.null(len))
      res = seq( dx[1], by=by, len=len)
    else if( !is.null(by))
      res= seq( dx[1], dx[2], by=by )
    else if( !is.null( len)) {
      if( length(x)> 1) 
         res = seq( dx[1], dx[2], len=len )
      else {
         res = seq( dx[1], dx[1]+len-1, by=1)
      }
    } else  
     res = seq( from=dx[1], to=dx[2], by=1 )
    sort(res)
}

rolling = function( a, dims=c(2), days=7 ) {
    if( is.null(dim(a))) {
        x = cumsum(a)
        c( rep(NA,days), tail(x,-days)-head(x,-days)) / days
    } else {
        apply(a, dims, function(x) {
        x = cumsum(x)
        c( rep(NA,days), tail(x,-days)-head(x,-days)) / days
        })
    }
}