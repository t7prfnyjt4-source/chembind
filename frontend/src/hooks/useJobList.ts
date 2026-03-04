import { useEffect, useMemo, useState } from "react";
import { collection, limit, onSnapshot, orderBy, query } from "firebase/firestore";
import { db } from "../firebase";
import type { JobDoc } from "./useBatchJob";

type UseJobListResult = {
  jobs: Array<JobDoc & { id: string }>;
  loading: boolean;
  error: string | null;
};

export function useJobList(uid: string | null | undefined): UseJobListResult {
  const [jobs, setJobs] = useState<Array<JobDoc & { id: string }>>([]);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);

  const enabled = useMemo(() => Boolean(uid), [uid]);

  useEffect(() => {
    if (!enabled) {
      setJobs([]);
      setLoading(false);
      setError(null);
      return;
    }

    setLoading(true);
    setError(null);

    const ref = collection(db, "users", uid!, "jobs");
    const q = query(ref, orderBy("createdAt", "desc"), limit(20));

    const unsub = onSnapshot(
      q,
      (snap) => {
        setJobs(snap.docs.map((d) => ({ id: d.id, ...(d.data() as JobDoc) })));
        setLoading(false);
      },
      (e) => {
        setError(e.message || "Failed to load job list.");
        setLoading(false);
      }
    );

    return () => unsub();
  }, [enabled, uid]);

  return { jobs, loading, error };
}
