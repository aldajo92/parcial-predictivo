function [REF] = ref_expand(ref,Hp)
    REF = kron(ones(Hp,1),ref);
end