function [check] = fun_check_adm(BARrows,DIAMrows,BARcolumns,DIAMcolumns,eta)
check=min([DIAMrows,DIAMcolumns])<eta*norm(BARrows-BARcolumns);
end

