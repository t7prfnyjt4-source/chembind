import type { AnalyzeRequest, AnalyzeResponse } from "./types";

const API_BASE = import.meta.env.VITE_API_BASE || "http://localhost:8000";

export class APIError extends Error {
  constructor(
    message: string,
    public status: number,
    public requestId?: string
  ) {
    super(message);
    this.name = "APIError";
  }
}

export async function analyzeSMILES(smiles: string): Promise<AnalyzeResponse> {
  const reqId = crypto.randomUUID();

  try {
    const res = await fetch(`${API_BASE}/analyze`, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ smiles } as AnalyzeRequest),
    });

    if (!res.ok) {
      throw new APIError(`Backend error: ${res.status}`, res.status, reqId);
    }

    return (await res.json()) as AnalyzeResponse;
  } catch (e) {
    if (e instanceof APIError) throw e;
    throw new APIError("Network error: Could not reach backend", 0, reqId);
  }
}
