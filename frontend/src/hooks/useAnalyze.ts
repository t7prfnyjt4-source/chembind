import { useCallback, useRef, useState } from "react";
import { analyzeSMILES, APIError } from "../api/client";
import type { AnalyzeResponse } from "../api/types";

export function useAnalyze() {
  const [data, setData] = useState<AnalyzeResponse | null>(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);

  // prevents race conditions: only latest request can update state
  const latestId = useRef<string | null>(null);

  const analyze = useCallback(async (smiles: string) => {
    const id = crypto.randomUUID();
    latestId.current = id;

    setLoading(true);
    setError(null);
    setData(null);

    try {
      const result = await analyzeSMILES(smiles);
      if (latestId.current === id) setData(result);
    } catch (err) {
      if (latestId.current === id) {
        setError(err instanceof APIError ? err.message : "Unknown error");
      }
    } finally {
      if (latestId.current === id) setLoading(false);
    }
  }, []);

  return { data, loading, error, analyze };
}
