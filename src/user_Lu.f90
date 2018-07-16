subroutine ReadHeatinject
include 'precision.inc'
include 'params.inc'

call AdvanceToNextInputLine( 4 ) ! add to the last line of the input file and param.inc                                                                         
read(4,*) iinj1, iinj2, jinj1, jinj2, xlatheat, ratfac

return
end
