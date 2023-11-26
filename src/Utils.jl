# Utils.jl
# 
# Contains general utility functions

################################################################################
# Type Conversion Utilities
################################################################################

# Time to # Seconds
secs = x -> Dates.value(Dates.Time(x)) / NANOSEC_TO_SEC
mins = x -> Dates.value(Dates.Time(x)) / NANOSEC_TO_SEC / SEC_TO_MIN

# Time to Index
time2index = (t,T) -> findfirst(x -> x == t, T)
index2time = (i,T) -> T[i]

mins2secs(x) = x * SEC_TO_MIN
secs2mins(x) = x / SEC_TO_MIN

################################################################################
# DataFrame Conversion Utilities
################################################################################

function dict2df(d, T, cols)
    df = DataFrame(time=T)
    for c in cols
        df[:,c] .= Real(0.0)
    end
    for (ti,t) in enumerate(T)
        for c in cols
            df[ti,c] = d[ti,c]
        end
    end
    return df
end
