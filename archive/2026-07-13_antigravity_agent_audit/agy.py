import asyncio
import os
import sys
import re
import math
from google.antigravity import Agent, LocalAgentConfig

def print_status(msg):
    """Helper function to print clean progress updates."""
    print(f"[*] {msg}")

async def execute_audit_with_retry(prompt, initial_model="gemini-3.5-flash"):
    """Executes the agent loop and handles rate limit recovery interactively."""
    current_model = initial_model
    
    while True:
        print_status(f"Preparing Agent using model: {current_model}")
        
        config = LocalAgentConfig(
            system_instructions="You are an expert Python codebase auditor. Be concise, act autonomously, and focus on actionable insights.",
            model=current_model 
        )
        
        try:
            async with Agent(config) as agent:
                print_status("Agent connected. Analyzing Git history and reading code...")
                
                response = await agent.chat(prompt)
                final_text = await response.text()
                
                print("\n" + "="*60)
                print("✅ FINAL AUDIT SUMMARY")
                print("="*60)
                print(final_text)
                
                return # Exit the loop on success
                
        except Exception as e:
            error_msg = str(e)
            
            # Catch rate limit and quota errors
            if "Quota exceeded" in error_msg or "429" in error_msg or "RESOURCE_EXHAUSTED" in error_msg:
                
                print("\n" + "!"*60)
                print("⚠️ RATE LIMIT REACHED")
                print("!"*60)
                
                # Print just the first sentence of the error so it remains readable
                clean_error = error_msg.split('\n')[0]
                print(f"Error detail: {clean_error}...") 
                
                # Present the generic options first
                print("\nChoose a recovery option:")
                print("  [1] Wait for the API's required delay time and retry")
                print("  [2] Downgrade to 'gemini-3.1-flash-lite' and retry immediately")
                print("  [3] Abort script")
                
                choice = input("\nEnter your choice (1/2/3): ").strip()
                
                if choice == '1':
                    # Automatically evaluate the regex and find the next integer delay
                    match = re.search(r"retry in\s*([\d\.]+)\s*s", error_msg, re.IGNORECASE)
                    
                    if match:
                        suggested_delay = math.ceil(float(match.group(1)))
                    else:
                        suggested_delay = 60 # Failsafe
                        
                    print_status(f"Waiting {suggested_delay}s...")
                    await asyncio.sleep(suggested_delay)
                    print_status("Resuming execution...")
                    
                elif choice == '2':
                    print_status("Downgrading model...")
                    current_model = "gemini-3.1-flash-lite"
                elif choice == '3':
                    print("Audit aborted by user.")
                    sys.exit(1)
                else:
                    print("Invalid choice. Aborting.")
                    sys.exit(1)
            else:
                # If it is a different, unrelated crash, print it and kill the script
                print(f"\n❌ A fatal error occurred: {e}")
                sys.exit(1)


async def run_clean_audit():
    print_status("Starting Antigravity audit...")
    
    if not os.environ.get("GEMINI_API_KEY"):
        print("❌ GEMINI_API_KEY environment variable is missing!")
        print("Please run 'export GEMINI_API_KEY=\"your_api_key\"' before running this script.")
        sys.exit(1)

    # Updated prompt to avoid shell pipes that trigger security blocks
    prompt = """
    Skip the Git commit history analysis. 
    
    Read and review the following 5 files for token optimization and algorithmic efficiency:
    1. compass.py
    2. damping_sweep_analysis.py
    3. avalanche_processor.py
    4. damping_sweep.py
    5. compass_generate_images.py
    
    CRITICAL INSTRUCTION: Do NOT delegate this task. Do NOT run this in the background. You must read these files and generate the structured findings into `audit_summary.json` IMMEDIATELY during this turn.
    
    Output a short, concise summary of your findings here once the file is written.
    """
    
    await execute_audit_with_retry(prompt)

if __name__ == "__main__":
    asyncio.run(run_clean_audit())
