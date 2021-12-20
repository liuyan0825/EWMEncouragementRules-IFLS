function f = predictp(X,Z,gamma)
% predict propensity scores based on logistic regression
XZ = [X Z];
f = exp(XZ*gamma)./(1+exp(XZ*gamma));
end