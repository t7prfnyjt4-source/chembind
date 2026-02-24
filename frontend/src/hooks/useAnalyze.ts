import { useState } from "react";
import { analyzeSMILES } from "../api/client";
import type { AnalyzeResponse } from "../api/types";

export function useAnalyze() {
  const [data, setData] = useState<AnalyzeResponse | null>(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);

  async function analyze(smiles: string) {
    setLoading(true);
    setError(null);
    setData(null);

    try {
      const result = await analyzeSMILES(smiles);
      setData(result);
      return result;
    } catch (e: any) {
      setError(e?.message ?? "Analyze failed");
      throw e;
    } finally {
      setLoading(false);
    }
  }

  return { data, loading, error, analyze };
}
