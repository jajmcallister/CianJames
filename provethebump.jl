# M will have a bump if M/I > m/i






function check_inequality(a, b, c, d)
    lhs = (a + b + c + d)^2 / 4
    rhs = a*d + c*d + b*c
    return lhs < rhs
end

using ProgressLogging

# Parameter sweep over ranges of a, b, c, d
@progress for a in 0.01:0.001:1
    for b in 0.01:0.001:1
        for c in 0.01:0.01:1
            for d in 0.01:0.01:1
                if check_inequality(a, b, c, d)
                    println("Inequality holds for a=$a, b=$b, c=$c, d=$d")
                end
            end
        end
    end
end