import { useCallback, useState } from "react";
import type { AnalyzeResponse } from "../api/types";
import { analyzeSMILES } from "../api/client";

type UseAnalyzeState = {
  loading: boolean;
  error: string | null;
  data: AnalyzeResponse | null;
};

export function useAnalyze() {
  const [state, setState] = useState<UseAnalyzeState>({
    loading: false,
    error: null,
    data: null,
  });

  const analyze = useCallback(async (smiles: string) => {
    setState((s) => ({ ...s, loading: true, error: null }));

    try {
      const res = await analyzeSMILES(smiles);
      setState({ loading: false, error: null, data: res });
      return res;
    } catch (e: any) {
      const msg =
        typeof e?.message === "string" && e.message.length > 0
          ? e.message
          : "Analyze failed";
      setState((s) => ({ ...s, loading: false, error: msg }));
      throw e;
    }
  }, []);

  return { ...state, analyze };
}
