from numpy import sum

def FSA(X,B):
    # assumes BOOZER coordinates
    return sum(X/B**2)/sum(1/B**2)
