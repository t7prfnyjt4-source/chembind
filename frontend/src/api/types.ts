export interface DescriptorResult {
  mw: number;
  logp: number;
  hbd: number;
  hba: number;
  tpsa: number;
  atom_count: number;
}

export interface AnalyzeResponse {
  apiVersion: string;
  requestId: string;
  smiles: string;
  canonicalSmiles: string | null;
  valid: boolean;
  descriptors: DescriptorResult | null;
  error: string | null;
}

export interface AnalyzeRequest {
  smiles: string;
}
