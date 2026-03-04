import { useEffect, useMemo, useState } from "react";
import { collection, limit, onSnapshot, orderBy, query } from "firebase/firestore";
import { db } from "../firebase";

export type JobItemDoc = {
  rowId: number;
  smiles: string;
  canonical?: string;
  valid?: boolean;

  mw?: number;
  logp?: number;
  hbd?: number;
  hba?: number;
  tpsa?: number;

  error?: string;
};

type UseJobItemsResult = {
  items: Array<JobItemDoc & { id: string }>;
  loading: boolean;
  error: string | null;
};

export function useJobItems(
  uid: string | null | undefined,
  jobId: string | null | undefined,
  enabledOverride?: boolean
): UseJobItemsResult {
  const [items, setItems] = useState<Array<JobItemDoc & { id: string }>>([]);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);

  const enabled = useMemo(
    () => Boolean(uid && jobId && (enabledOverride ?? true)),
    [uid, jobId, enabledOverride]
  );

  useEffect(() => {
    if (!enabled) {
      setItems([]);
      setLoading(false);
      setError(null);
      return;
    }

    setLoading(true);
    setError(null);

    const ref = collection(db, "users", uid!, "jobs", jobId!, "items");
    const q = query(ref, orderBy("rowId", "asc"), limit(2000));

    const unsub = onSnapshot(
      q,
      (snap) => {
        setItems(snap.docs.map((d) => ({ id: d.id, ...(d.data() as JobItemDoc) })));
        setLoading(false);
      },
      (e) => {
        setError(e.message || "Failed to load job items.");
        setLoading(false);
      }
    );

    return () => unsub();
  }, [enabled, uid, jobId]);

  return { items, loading, error };
}
