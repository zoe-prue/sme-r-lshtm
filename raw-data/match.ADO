*! version 2.0.0 mh 9 june 2002
program define match
version 8.0
syntax varlist(min=3 max=3) [if] [in]

preserve

if "`if'"!="" | "`in'"!="" {
        qui keep `if' `in'
}

collapse (sum) "`2'" , by("`1'" "`3'") cw

qui reshape wide `2', i(`3') j(`1')
qui gen _cont=`2'0
qui gen _case=`2'1
di in gr _n(1) "Exposure status of case and controls in the same set"
di ""
tab _cont _case



restore 
end

