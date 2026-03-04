import { useEffect, useMemo, useState } from "react";
import { doc, onSnapshot } from "firebase/firestore";
import { db } from "../firebase"; // adjust if your firebase file path differs

export type JobStatus = "queued" | "running" | "completed" | "failed";

export type JobDoc = {
  status: JobStatus;
  createdAt?: any;
  total?: number;
  processed?: number;
  successCount?: number;
  failureCount?: number;
  errorMessage?: string;
};

type UseBatchJobResult = {
  job: (JobDoc & { id: string }) | null;
  loading: boolean;
  error: string | null;
};

export function useBatchJob(
  uid: string | null | undefined,
  jobId: string | null | undefined
): UseBatchJobResult {
  const [job, setJob] = useState<(JobDoc & { id: string }) | null>(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);

  const enabled = useMemo(() => Boolean(uid && jobId), [uid, jobId]);

  useEffect(() => {
    if (!enabled) {
      setJob(null);
      setLoading(false);
      setError(null);
      return;
    }

    setLoading(true);
    setError(null);

    const ref = doc(db, "users", uid!, "jobs", jobId!);

    const unsub = onSnapshot(
      ref,
      (snap) => {
        if (!snap.exists()) {
          setJob(null);
          setError("Job not found.");
          setLoading(false);
          return;
        }
        setJob({ id: snap.id, ...(snap.data() as JobDoc) });
        setLoading(false);
      },
      (e) => {
        setError(e.message || "Failed to subscribe to job.");
        setLoading(false);
      }
    );

    return () => unsub();
  }, [enabled, uid, jobId]);

  return { job, loading, error };
}
