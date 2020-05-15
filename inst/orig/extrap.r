extrap=function(efac,x,x1,x2,y1,y2) {
   if(efac<0){ x=pmax(x,x1+efac*(x2-x1))
      y1+(y2-y1)*(x-x1)/(x2-x1)
   }
   else { x=pmin(x,x2+efac*(x2-x1))
      y2+(y2-y1)*(x-x2)/(x2-x1)}
}
