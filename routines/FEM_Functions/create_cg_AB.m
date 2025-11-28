function [A,B]=create_cg_AB(mesh,func_data)
    
    Txx=func_data.Txx; Tyy=func_data.Tyy; B=func_data.T;
    
    A=Txx+Tyy;
end