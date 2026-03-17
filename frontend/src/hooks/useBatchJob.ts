import { useEffect, useMemo, useState } from "react";
import { doc, onSnapshot, Timestamp } from "firebase/firestore";
import { db } from "../firebase/config";

export type JobStatus = "queued" | "running" | "completed" | "failed";

export type JobDoc = {
  status: JobStatus;
  createdAt?: Timestamp | { seconds: number } | null;
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

function errMsg(e: unknown): string {
  if (e instanceof Error) return e.message;
  return String(e);
}

export function useBatchJob(
  uid: string | null | undefined,
  jobId: string | null | undefined
): UseBatchJobResult {
  const [job, setJob] = useState<(JobDoc & { id: string }) | null>(null);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);

  const enabled = useMemo(() => Boolean(uid && jobId), [uid, jobId]);

  useEffect(() => {
    if (!enabled) return;

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
        setError(null);
        setLoading(false);
      },
      (e: unknown) => {
        setError(errMsg(e) || "Failed to subscribe to job.");
        setLoading(false);
      }
    );

    return () => unsub();
  }, [enabled, uid, jobId]);

  return { job, loading, error };
}
