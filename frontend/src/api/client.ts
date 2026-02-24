import type { AnalyzeRequest, AnalyzeResponse } from "./types";
import { auth } from "../firebase/config";

function getApiBase(): string {
  const raw = import.meta.env.VITE_API_BASE || "http://localhost:8000";
  // Trim trailing slashes to avoid // in URLs
  return String(raw).replace(/\/+$/, "");
}

export class APIError extends Error {
  status: number;
  requestId: string;

  constructor(message: string, status: number, requestId: string) {
    super(message);
    this.name = "APIError";
    this.status = status;
    this.requestId = requestId;
  }
}


type RequestOptions = {
  method?: "GET" | "POST" | "PUT" | "DELETE";
  headers?: Record<string, string>;
  body?: unknown;
  idempotencyKey?: string;

  // fail-safe knobs
  timeoutMs?: number; // default 20000
  retries?: number; // default 1
  retryDelayMs?: number; // default 400
};

async function sleep(ms: number) {
  await new Promise((r) => setTimeout(r, ms));
}

async function readErrorBody(res: Response): Promise<string> {
  try {
    const ct = res.headers.get("content-type") || "";
    if (ct.includes("application/json")) return JSON.stringify(await res.json());
    return await res.text();
  } catch {
    return "";
  }
}

export async function apiFetch<T>(path: string, opts: RequestOptions = {}): Promise<T> {
  const API_BASE = getApiBase();
  const reqId = crypto.randomUUID();

  const timeoutMs = opts.timeoutMs ?? 20_000;
  const retries = opts.retries ?? 1;
  const retryDelayMs = opts.retryDelayMs ?? 400;

  const url = `${API_BASE}${path.startsWith("/") ? "" : "/"}${path}`;

  const user = auth.currentUser;
  const token = user ? await user.getIdToken() : null;

  const headers: Record<string, string> = {
    "Content-Type": "application/json",
    ...(opts.headers ?? {}),
  };
  if (token) headers["Authorization"] = `Bearer ${token}`;
  if (opts.idempotencyKey) headers["Idempotency-Key"] = opts.idempotencyKey;

  let attempt = 0;
  while (true) {
    attempt += 1;

    const controller = new AbortController();
    const t = setTimeout(() => controller.abort(), timeoutMs);

    try {
      const res = await fetch(url, {
        method: opts.method ?? "GET",
        headers,
        body: opts.body !== undefined ? JSON.stringify(opts.body) : undefined,
        signal: controller.signal,
      });

      if (!res.ok) {
        const msg = await readErrorBody(res);
        const err = new APIError(
          `Backend error: ${res.status}${msg ? `: ${msg}` : ""}`,
          res.status,
          reqId
        );

        // Retry only on 5xx
        if (res.status >= 500 && attempt <= retries) {
          await sleep(retryDelayMs * attempt);
          continue;
        }
        throw err;
      }

      const ct = res.headers.get("content-type") || "";
      if (!ct.includes("application/json")) {
        return (undefined as unknown) as T;
      }
      return (await res.json()) as T;
    } catch (e: any) {
      const isAbort = e?.name === "AbortError";
      const isNetwork = !isAbort && (e instanceof TypeError || String(e?.message || "").includes("Failed to fetch"));

      if ((isAbort || isNetwork) && attempt <= retries) {
        await sleep(retryDelayMs * attempt);
        continue;
      }

      if (e instanceof APIError) throw e;
      throw new APIError(
        isAbort ? `API timeout after ${timeoutMs}ms` : "Network error: Could not reach backend",
        0,
        reqId
      );
    } finally {
      clearTimeout(t);
    }
  }
}

export async function analyzeSMILES(smiles: string): Promise<AnalyzeResponse> {
  return apiFetch<AnalyzeResponse>("/api/analyze", {
    method: "POST",
    body: ({ smiles } as AnalyzeRequest),
    retries: 1,
    timeoutMs: 20_000,
  });
}
