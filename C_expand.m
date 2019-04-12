function [C_exp] = C_expand(C,Hp)
    I_C = eye(Hp);
    C_exp = kron(I_C,C);
end