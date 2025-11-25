function output = I_set_interval_precision( var, precision )
   global INTERVAL_MODE;
   if INTERVAL_MODE
    output = var + infsup(-precision, precision);
   else
    output = var;
   end
   INTERVAL_MODE
end
