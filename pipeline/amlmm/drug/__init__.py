"""COMPASS-AML — the ex-vivo drug-response layer of MOSAIC-AML.

Predicts, for an uploaded patient, which small-molecule inhibitors their leukaemia is likely to be
sensitive to, learned from the BeatAML2 ex-vivo functional screen. Deliberately kept parallel to (not
downstream of) the mutation caller: mutation calls enter as *evidence*, never as a drug look-up.

  data        BeatAML2 probit curve fits -> QC'd, within-drug normalised response table
  targets     curated inhibitor -> target / mechanism / clinical-tier annotation
  features    patient feature blocks (bulk RNA, mutations, clinical, differentiation state)
  model       Model A - hierarchical patient-level response model (regression + tail classifier)
  statemodel  Model B - the same model applied per malignant cell state -> coverage metrics
  mechanism   Model C - target-pathway activity, independent of the empirical response model
  utility     the treatment-utility score S_ij and its tiered rankings
  agents      the eight drug-reasoning expert agents
"""
