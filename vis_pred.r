for (i in 1:N) {
  plot(1:(Tall+7),uk_cases[i,3:(Tall+2+7)],main=uk_cases[i,2])
  lines(Tcond+Tlik+(1:Tproj),Cprojlower[i,])
  lines(Tcond+Tlik+(1:Tproj),Cprojmedian[i,])
  lines(Tcond+Tlik+(1:Tproj),Cprojupper[i,])
  # b <- scan("stdin",character(), n=1)
  invisible(readline())
}
  
