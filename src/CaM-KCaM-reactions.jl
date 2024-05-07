"""
$(SIGNATURES)

The reactions described in this function were taken from Chang et al. 2019.
The tag simpler refers to the simplest implementation proposed by Chang.
Note also that Calcineurin is added to reactions below.
"""
ChangCaMKII_simpler = @reaction_network begin
  @parameters kf_2C kb_2C kf_2N kb_2N kf_CaM0 kb_CaM0 kf_CaM2N kb_CaM2N kf_CaM2C kb_CaM2C kf_CaM4 kb_CaM4 kf_K2N kb_K2N kf_K2C kb_K2C F k1 k2 k4 k5 k3 kcanf kcanb
  (kf_2C,kb_2C), CaM0 + 2Ca ↔ CaM2C
  (kf_2N,kb_2N), CaM0 +  2Ca ↔ CaM2N
  (kf_2N,kb_2N), CaM2C + 2Ca ↔ CaM4
  (kf_2C,kb_2C), CaM2N + 2Ca ↔ CaM4

  (kcanf,kcanb), CaM4 + mCaN ↔ CaN4

  (kf_CaM0,kb_CaM0), CaM0 + mKCaM ↔ KCaM0
  (kf_CaM2N,kb_CaM2N), CaM2N + mKCaM ↔ KCaM2N
  (kf_CaM2C,kb_CaM2C), CaM2C + mKCaM ↔ KCaM2C
  (kf_CaM4,kb_CaM4), CaM4 + mKCaM ↔ KCaM4

  (kf_K2N,kb_K2N), KCaM0 + 2Ca ↔ KCaM2N
  (kf_K2C,kb_K2C), KCaM0 + 2Ca ↔ KCaM2C
  (kf_K2N,kb_K2N), KCaM2C + 2Ca ↔ KCaM4
  (kf_K2C,kb_K2C), KCaM2N + 2Ca ↔ KCaM4

  (F*k1), KCaM0 →  PCaM0
  (F*k1), KCaM2N →  PCaM2N
  (F*k1), KCaM2C →  PCaM2C
  (F*k1), KCaM4 →  PCaM4

  (k2), PCaM0 → P + CaM0
  (k2), PCaM2N → P + CaM2N
  (k2), PCaM2C → P + CaM2C
  (k2), PCaM4 → P + CaM4

  (k4,k5), P ↔ P2
  (k3), P → mKCaM

end

"""
$(SIGNATURES)

This function extracts the differential equations from the CaM-CaMKII-CaN reactions and write it into the file with name `filename = "write-equation.txt"`. The outcome of this function is used to write the `F_synapse`.
"""
function extractReactionsCamKII(filename = "write-equation.txt")
  f = open(filename, "w")
  rn = ChangCaMKII_simpler
  osys = convert(ODESystem, rn)
  eqns = equations(osys)
  for (i, eq) in enumerate(eqns)
      write(f, @sprintf("%s=%s \n", i, eq))
  end
  close(f)
end

function extractReactionsCamKII_prob()
  u0_symbol = [:CaM0 => 1., :Ca => 1., :CaM2C => 1., :CaM2N => 1., :CaM4 => 1., :mCaN => 1., :CaN4 => 1., :mKCaM => 1., :KCaM0 => 1., :KCaM2N => 1., :KCaM2C => 1., :KCaM4 => 1., :PCaM0 => 1., :PCaM2N => 1., :PCaM2C => 1., :PCaM4 => 1., :P => 1., :P2 => 1., ]

  u0 = rand(length(u0_symbol))

  par = (:kf_2C => 1.0, :kb_2C => 1.0, :kf_2N => 1.0, :kb_2N => 1.0, :kf_CaM0 => 1.0, :kb_CaM0 => 1.0, :kf_CaM2N => 1.0, :kb_CaM2N => 1.0, :kf_CaM2C => 1.0, :kb_CaM2C => 1.0, :kf_CaM4 => 1.0, :kb_CaM4 => 1.0, :kf_K2N => 1.0, :kb_K2N => 1.0, :kf_K2C => 1.0, :kb_K2C => 1.0, :F => 1.0, :k1 => 1.0, :k2 => 1.0, :k4 => 1.0, :k5 => 1.0, :k3 => 1.0, :kcanf => 1.0, :kcanb => 1.0)

  oprob = ODEProblem(ChangCaMKII_simpler, u0, (0., 1.), par)
end
