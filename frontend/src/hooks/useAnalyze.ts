import { useState } from "react";
import { auth } from "../firebase.ts";


const API_BASE = import.meta.env.VITE_API_BASE;

export function useAnalyze() {
  const [data, setData] = useState<any>(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);

  async function analyze(smiles: string) {
    setLoading(true);
    setError(null);
    setData(null);

    try {
      const user = auth.currentUser;
      if (!user) {
        throw new Error("Not authenticated");
      }

      const token = await user.getIdToken();

      const res = await fetch(`${API_BASE}/analyze`, {
        method: "POST",
        headers: {
          "Content-Type": "application/json",
          "Authorization": `Bearer ${token}`,
        },
        body: JSON.stringify({ smiles }),
      });

      if (!res.ok) {
        const text = await res.text();
        throw new Error(text);
      }

      const json = await res.json();
      setData(json);
    } catch (e: any) {
      setError(e.message || "Unknown error");
    } finally {
      setLoading(false);
    }
  }

  return { data, loading, error, analyze };
}



